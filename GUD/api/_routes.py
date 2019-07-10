from GUD.api import app
from GUD.api.db import establish_GUD_session, shutdown_session
from flask import request, jsonify
from GUD.ORM import (Gene, ShortTandemRepeat, CNV, ClinVar, Conservation,
                     DNAAccessibility, Enhancer, HistoneModification, TAD, TFBinding, TSS)
from GUD.ORM.genomic_feature import GenomicFeature
from werkzeug.exceptions import HTTPException, NotFound, BadRequest
import sys
import math
import re
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
        uids = uids.split(',')
        for i in range(len(uids)):
            if uids[i].isdigit():
                uids[i] = int(uids[i])
        if len(uids) > 1000 or len(uids) < 1:
            raise BadRequest(
                "list of uids must be greater than 0 and less than 1000")
        return resource.select_by_uids(session, uids, limit, offset)
    elif {'chrom', 'end', 'location', 'sources', 'start'} == keys or {'chrom', 'end', 'location', 'start'} == keys:
        try:
            start = int(start) - 1
            end = int(end)
        except:
            raise BadRequest("start and end should be formatted as integers")
        if re.fullmatch('^chr(X|Y|[1-9]|1[0-9]|2[0-2])$', chrom) == None:
            raise BadRequest("chromosome should be formatted as chrZ where z \
                is X, Y, or 1-22")
        if (end - start) > 4000000:
            raise BadRequest("start and end must be less than 4000000bp apart")
        if location not in ['within', 'overlapping', 'exact']:
            raise BadRequest("location must be specified as withing, overlapping, or exact")
        if 'sources' in keys:
            sources = request.args.get('sources', default=None)
            sources = sources.split(',')
            return resource.select_by_sources(session, chrom, start, end, location, sources, limit, offset)
        return resource.select_by_location(session, chrom, start, end, location, limit, offset)


def genomic_feature_mixin2_queries(session, resource, request, limit, offset):
    keys = set(request.args)
    keys.discard('page')
    chrom = request.args.get('chrom', default=None, type=str)
    end = request.args.get('end', default=None)
    start = request.args.get('start', default=None)
    location = request.args.get('location', default=None)

    if {'chrom', 'end', 'location', 'samples', 'start'} == keys or \
            {'chrom', 'end', 'experiments', 'location', 'start'} == keys:
        try:
            start = int(start) - 1
            end = int(end)
        except:
            raise BadRequest("start and end should be formatted as integers")
        if re.fullmatch('^chr(X|Y|[1-9]|1[0-9]|2[0-2])$', chrom) == None:
            raise BadRequest("chromosome should be formatted as chrZ where z \
                is X, Y, or 1-22")
        if (end - start) > 4000000:
            raise BadRequest("start and end must be less than 4000000bp apart")
        if location not in ['within', 'overlapping', 'exact']:
            raise BadRequest("location must be specified as withing, overlapping, or exact")
        if 'samples' in keys:
            samples = request.args.get('samples', default=None)
            samples = samples.split(',')
            samples = [s.replace("+", " ") for s in samples]
            if len(samples) > 1000 or len(samples) < 1:
                raise BadRequest(
                    "list of samples must be greater than 0 and less than 1000")
            return resource.select_by_samples(session, chrom, start, end, samples, location,  limit, offset)
        elif 'experiments' in keys:
            experiments = request.args.get('experiments', default=None)
            experiments = experiments.split(',')

            if len(experiments) > 1000 or len(experiments) < 1:
                raise BadRequest(
                    "list of samples must be greater than 0 and less than 1000")
            return resource.select_by_experiments(session, chrom, start, end, experiments, location,  limit, offset)
    else:
        return genomic_feature_mixin1_queries(session, resource, request, limit, offset)


def gene_queries(session, resource, request, limit, offset):
    keys = set(request.args)
    keys.discard('page')
    if {'names'} == keys:
        names = request.args.get('names', default=None)
        names = names.split(',')
        if len(names) > 1000 or len(names) < 1:
            raise BadRequest(
                "list of names must be greater than 0 and less than 1000")
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
    if {'motif'} == keys or {'motif', 'rotation'} == keys:
        motif = request.args.get('motif', default=None)
        rotation = request.args.get('rotation', default=False, type=bool)
        return resource.select_by_motif(session, motif.upper(), limit, offset, rotation)
    else:
        return genomic_feature_mixin1_queries(session, resource, request, limit, offset)


def tad_queries(session, resource, request, limit, offset):
    keys = set(request.args)
    keys.discard('page')
    if {'restriction_enzymes'} == keys:
        restriction_enzymes = request.args.get('restriction_enzymes', default=None)
        restriction_enzymes = restriction_enzymes.split(',')
        if len(restriction_enzymes) > 1000 or len(restriction_enzymes) < 1:
            raise BadRequest(
                "list of restriction enzymes must be greater than 0 and less than 1000")
        return resource.select_by_restriction_enzymes(session, restriction_enzymes, limit, offset)
    else:
        return genomic_feature_mixin2_queries(session, resource, request, limit, offset)


def tf_binding_queries(session, resource, request, limit, offset):
    keys = set(request.args)
    keys.discard('page')
    if {'restriction_enzymes'} == keys:
        tf = request.args.get('restriction_enzymes', default=None)
        return resource.select_by_tf(cls, session, chrom, start, end, location, tf, limit, offset)
    else:
        return genomic_feature_mixin2_queries(session, resource, request, limit, offset)


def tss_queries(session, resource, request, limit, offset):
    return false


@app.route('/api/v1/genes/symbols')
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

@app.route('/api/v1/short_tandem_repeats/pathogenic')
def pathogenic_str():
    page = request.args.get('page', default=1, type=int)
    if (page <= 0):
        raise BadRequest('pages must be positive integers')
    session     = establish_GUD_session()
    offset      = (page-1)*20
    limit       = 20
    resource    = ShortTandemRepeat()
    result      =  resource.select_by_pathogenicity(session, limit, offset)
    shutdown_session(session)
    result      = create_page(resource, result, page, request.url)
    return jsonify(result)


def check_page(request):
    page = request.args.get('page', default=1, type=int)
    if (page <= 0):
        raise BadRequest('pages must be positive integers')


@app.route('/api/v1/dna_accessibility')
def dna_accessibility():
    check_page(request)
    session = establish_GUD_session()
    dna_accessibility = DNAAccessibility()

    shutdown_session(session)
    result = create_page(resource, result, page, request.url)
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
    elif resource in ['dna_accessibility', 'histone_modifications', 'enhancers']:
        if resource == 'dna_accessibility':
            resource = DNAAccessibility()
        elif resource == 'histone_modifications':
            resource = HistoneModification()
        elif resource == 'enhancers':
            resource = Enhancer()
        result = genomic_feature_mixin2_queries(
            session, resource, request, limit, offset)
    elif resource == 'tads':
        resource = TAD()
        result = tad_queries(session, resource, request, limit, offset)
    elif resource == 'tf_binding':
        resource = TFBinding()
        result = tf_binding_queries(session, resource, request, limit, offset)
    elif resource == 'tss':
        resource = TSS()
        result = tss_queries(session, resource, request, limit, offset)
    else:
        raise BadRequest('valid resources are genes, short_tandem_repeats,\
             copy_number_variants, clinvar, conservation')

    shutdown_session(session)
    result = create_page(resource, result, page, request.url)
    return jsonify(result)
