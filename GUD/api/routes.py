from GUD.api import app
from GUD.api.db import establish_GUD_session, shutdown_session
from flask import request, jsonify
from GUD.ORM import Gene, ShortTandemRepeat, CNV, ClinVar, Conservation, DNAAccessibility
from GUD.ORM.genomic_feature import GenomicFeature
from werkzeug.exceptions import HTTPException, NotFound, BadRequest
import sys, math, re
# print(names, file=sys.stdout)


@app.route('/')
def index():
    return 'HOME'

def create_page(resource, result, page, url) -> dict:
    """
    returns 404 error or a page
    """

    page_size = 20
    result_size = result[0]
    json = {}
    if result_size == 0:
        raise NotFound('No results from this query')
    if page <= 0 or (page-1)*page_size > result_size:
        raise BadRequest('Page range is invalid, valid range for this query is between 1 and ' +
                         str(math.ceil(result_size/page_size)-1))

    results = result[1]

    if (resource != False):
        results = result[1]
        results = [resource.as_genomic_feature(e) for e in results]
        results = [e.serialize() for e in results]

    json = {'size': result_size,
            'results': results}
    if (page)*page_size < result_size:  # has next
        if re.search('\?', url) is None:
            next_page = url+'?page='+str(page+1)
        elif re.search('page', url) is None:
            next_page = url+'&page='+str(page+1)
        else:
            next_page = re.sub('page=\d+', 'page='+str(page+1), url)
        json['next'] = next_page
    if (page-2)*page_size >= 0:  # has prev
        prev_page = re.sub('page=\d+', 'page='+str(page-1), url)
        json['prev'] = prev_page
    return json


def genomic_feature_mixin1_queries(session, resource, request, limit, offset):
    keys = set(request.args)
    keys.discard('page')

    chrom = request.args.get('chrom', default=None, type=str)
    end = request.args.get('end', default=None)
    location = request.args.get('location', default=None, type=str)
    start = request.args.get('start', default=None)

    if {'uids'} == keys:
        uids = request.args.get('uids', default=None)
        try:
            uids = uids.split(',')
            uids = [int(e) for e in uids]
        except:
            raise BadRequest(
                "uids must be positive integers seperated by commas (,).")
        if len(uids) > 1000 or len(uids) < 1:
             raise BadRequest("list of uids must be greater than 0 and less than 1000")
        result = resource.select_by_uids(session, uids, limit, offset)
    elif {'chrom', 'end', 'sources', 'start'} == keys:
        sources = request.args.get('sources', default=None)
        sources = sources.split(',')
        try:
            start = int(start) - 1
            end = int(end)
        except:
            raise BadRequest("start and end should be formatted as integers, \
            chromosomes should be formatted as chrZ.")
        if (end - start) > 1000000:
            raise BadRequest("start and end must be less than 1000000bp apart")
        if len(sources) > 1000 or len(sources) < 1:
             raise BadRequest("list of sources must be greater than 0 and less than 1000")       
        result = resource.select_by_sources(session, chrom, start, end, sources, limit, offset)
    elif {'chrom', 'end', 'location', 'start'} == keys:
        try:
            start = int(start) - 1
            end = int(end)
        except:
            raise BadRequest("start and end should be formatted as integers, \
            chromosomes should be formatted as chrZ.")
        if (end - start) > 1000000:
            raise BadRequest("start and end must be less than 1000000bp apart")
        if location == 'exact':
            result = resource.select_by_exact_location(
                session, chrom, start, end, limit, offset)
        else:
            result = resource.select_by_location(
                session, chrom, start, end, limit, offset)
    else:
        raise BadRequest(
            "invalid combination of parameters for this resource, refer \
                    documentation for correct combination of paramaters")

    return result


def genomic_feature_mixin2_queries(session, resource, request, limit, offset):
    keys = set(request.args)
    keys.discard('page')
    chrom = request.args.get('chrom', default=None, type=str)
    end = request.args.get('end', default=None)
    start = request.args.get('start', default=None)
    
    if {'chrom', 'end', 'samples', 'start'} == keys:
        samples = request.args.get('samples', default=None)
        samples = samples.split(',')
        samples = [s.replace("+", " ") for s in samples]
        try:
            start = int(start) - 1
            end = int(end)
        except:
            raise BadRequest("start and end should be formatted as integers, \
            chromosomes should be formatted as chrZ.")
        if (end - start) > 1000000:
            raise BadRequest("start and end must be less than 1000000bp apart")
        if len(samples) > 1000 or len(samples) < 1:
             raise BadRequest("list of samples must be greater than 0 and less than 1000")       
        return resource.select_by_samples(session, chrom, start, end, samples, limit, offset)
    elif {'chrom', 'end', 'experiments', 'start'} == keys:
        experiments = request.args.get('experiments', default=None)
        experiments = experiments.split(',')
        try:
            start = int(start) - 1
            end = int(end)
        except:
            raise BadRequest("start and end should be formatted as integers, \
            chromosomes should be formatted as chrZ.")
        if (end - start) > 1000000:
            raise BadRequest("start and end must be less than 1000000bp apart")
        if len(experiments) > 1000 or len(experiments) < 1:
             raise BadRequest("list of samples must be greater than 0 and less than 1000") 
        return resource.select_by_experiments(session, chrom, start, end, experiments, limit, offset)
    
    else: 
        return genomic_feature_mixin1_queries(session, resource, request, limit, offset)


def gene_queries(session, resource, request, limit, offset):
    keys = set(request.args)
    keys.discard('page')
    if {'names'} == keys:
        names = request.args.get('names', default=None)
        names = names.split(',')
        if len(names) > 1000 or len(names) < 1:
             raise BadRequest("list of names must be greater than 0 and less than 1000")
        return resource.select_by_names(session, limit, offset, names)
    else:
        return genomic_feature_mixin1_queries(session, resource, request, limit, offset)


def clinvar_queries(session, resource, request, limit, offset):
    keys = set(request.args)
    keys.discard('page')
    if {'clinvar_id'} == keys:
        clinvarID = request.args.get('clinvar_id')
        return resource.select_by_clinvarID(session, clinvarID, limit, offset)
    else:
        return genomic_feature_mixin1_queries(session, resource, request, limit, offset)


def short_tandem_repeat_queries(session, resource, request, limit, offset):
    keys = set(request.args)
    keys.discard('page')
    if {'pathogenicity'} == keys:
        pathogenicity = request.args.get(
            'pathogenicity', default=None, type=bool)
        return resource.select_by_pathogenicity(session, limit, offset)
    elif {'motif'} == keys or {'motif', 'rotation'} == keys:
        motif = request.args.get('motif', default=None)
        rotation = request.args.get('rotation', default=False, type=bool)
        return resource.select_by_motif(session, motif.upper(), limit, offset, rotation)
    else:
        return genomic_feature_mixin1_queries(session, resource, request, limit, offset)


@app.route('/api/v1/genesymbols')
def gene_symbols():
    url = request.url
    session = establish_GUD_session()
    page = int(request.args.get('page', default=1))
    offset = (page-1)*20
    limit = 20
    result = Gene().get_all_gene_symbols(session, limit, offset)
    shutdown_session(session)
    result = create_page(False, result, page, url)
    return jsonify(result)


@app.route('/api/v1/<resource>')
def resource(resource):
    page = request.args.get('page', default=1, type=int)
    if (page <= 0):
        raise BadRequest('pages must be positive integers')
    session = establish_GUD_session()
    offset = (page-1)*20
    limit = 20
    result = None

    if resource == 'clinvar':
        resource = ClinVar()
        result = clinvar_queries(session, resource, request, limit, offset)
    elif resource == 'concervation':
        resource = Conservation()
        result = genomic_feature_mixin1_queries(
            session, resource, request, limit, offset)
    elif resource == 'copy_number_variants':
        resource = CNV()
        result = genomic_feature_mixin1_queries(
            session, resource, request, limit, offset)
    elif resource == 'genes':
        resource = Gene()
        result = gene_queries(session, resource, request, limit, offset)
    elif resource == 'short_tandem_repeats':
        resource = ShortTandemRepeat()
        result = short_tandem_repeat_queries(
            session, resource, request, limit, offset)
    elif resource == 'dna_accessibility':
        resource = DNAAccessibility()
        result = genomic_feature_mixin2_queries(
                session, resource, request, limit, offset)
    else:
        raise BadRequest('valid resources are genes, short_tandem_repeats,\
             copy_number_variants, clinvar, conservation')

    shutdown_session(session)
    result = create_page(resource, result, page, request.url)
    return jsonify(result)


# examples
# genes
# http://127.0.0.1:5000/api/v1/genesymbols
# http://127.0.0.1:5000/api/v1/genes?uids=1
# http://127.0.0.1:5000/api/v1/genes?names=LOC102725121
# http://127.0.0.1:5000/api/v1/genes?chrom=chr1&start=11869&end=14362&location=exact ######
# http://127.0.0.1:5000/api/v1/genes?sources=refGene&chrom=chr1&start=11869&end=14362
# str
# http://127.0.0.1:5000/api/v1/short_tandem_repeats?pathogenicity=True
# http://127.0.0.1:5000/api/v1/short_tandem_repeats?motif=ATGGG&rotation=True
# CNV
# http://127.0.0.1:5000/api/v1/copy_number_variants?chrom=chr1&start=120672583&end=142552733&location=exact
# Clinvar
# http://127.0.0.1:5000/api/v1/clinvar?clinvar_id=475283
# http://127.0.0.1:5000/api/v1/clinvar?uids=1
# http://127.0.0.1:5000/api/v1/clinvar?chrom=chr1&start=949422&end=949422&location=exact
# http://127.0.0.1:5000/api/v1/clinvar?sources=ClinVar_2018-10-28&chrom=chr1&start=11869&end=14362
# dna_accessibility
# http://127.0.0.1:5000/api/v1/dna_accessibility?uids=1
# http://127.0.0.1:5000/api/v1/dna_accessibility?chrom=chr1&start=10410&end=10606&location=exact 
# http://127.0.0.1:5000/api/v1/dna_accessibility?chrom=chr1&start=1&end=100000&location=within
# http://127.0.0.1:5000/api/v1/dna_accessibility?sources=ENCODE&chrom=chr1&start=1&end=100000
# http://127.0.0.1:5000/api/v1/dna_accessibility?samples=mesoderm+(heart)&chrom=chr1&start=11869&end=14362
# http://127.0.0.1:5000/api/v1/dna_accessibility?chrom=chr1&start=1&end=100000&experiments=DNase-seq