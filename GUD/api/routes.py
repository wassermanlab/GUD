from GUD.api import app
from GUD.api.db import establish_GUD_session, shutdown_session
from flask import request, jsonify
from GUD.ORM import Gene, ShortTandemRepeat
from GUD.ORM.genomic_feature import GenomicFeature
import re
from werkzeug.exceptions import HTTPException, NotFound, BadRequest
import sys, math
# print(names, file=sys.stdout)
# errors
# 400 bad request, client input validation fails
# 404 not found


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
        raise BadRequest('Page range is invalid, valid range for this query is between 1 and ' + str(math.ceil(result_size/page_size)-1))
    
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


def genomic_feature_queries(session, resource, uids, chrom, start, end, sources, location, limit, offset): 
    if uids is not None and all(v is None for v in [chrom, start, end, sources, location]):
        try:
            uids = uids.split(',')
            uids = [int(e) for e in uids]
        except:
            raise BadRequest(
                "uids must be positive integers seperated by commas (,).")
        result = resource.select_by_uids(session, uids, limit, offset)
    elif sources is not None and all(v is None for v in [uids, chrom, start, end, location]):
        sources = sources.split(',')
        result = resource.select_by_sources(session, sources, limit, offset)
    elif chrom is not None and start is not None and end is not None and all(v is None for v in [uids, sources]):
        try:
            start = int(start) - 1
            end = int(end)
        except:
            raise BadRequest("start and end should be formatted as integers, \
            chromosomes should be formatted as chrZ.")
        if location == 'exact':
            result = resource.select_by_exact_location(
                session, chrom, start, end, limit, offset)
        else:
            result = resource.select_by_location(session, chrom, start, end, limit, offset)
    else: 
        raise BadRequest('requests must have some parameters, refer to the \
            docs for the correct parameters')

    return result


def gene_queries(session, resource, names, limit, offset):
    names = names.split(',')
    return resource.select_by_names(session, limit, offset, names)


def short_tandem_repeat_queries(session, resource, pathogenicity, motif, rotation, limit, offset):
    if pathogenicity is True and all(v is None for v in [motif, rotation]):
        print(pathogenicity, file=sys.stdout)
        return resource.select_by_pathogenicity(session, limit, offset)
    elif rotation is True:
        return resource.select_by_motif(session, motif.upper(), limit, offset, rotation)
    else: 
        return resource.select_by_motif(session, motif.upper(), limit, offset)


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
    session = establish_GUD_session()
    # parameters
    page = request.args.get('page', default=1, type=int)
    if (page <= 0):
        raise BadRequest('pages must be positive integers')
    offset = (page-1)*20
    limit = 20
    uids = request.args.get('uids', default=None)
    chrom = request.args.get('chrom', default=None)
    start = request.args.get('start', default=None)
    end = request.args.get('end', default=None)
    sources = request.args.get('sources', default=None)
    location = request.args.get('location', default=None, type=str)
    result = None

    #queries unique to resources
    if (resource == 'genes'):
        names = request.args.get('names', default=None)
        resource = Gene()
        if names is not None and all(v is None for v in [uids, chrom, start, end, sources, location]):
            result = gene_queries(session, resource, names, limit, offset)

    elif (resource == 'short_tandem_repeats'):
        pathogenicity   = request.args.get('pathogenicity', default=None, type=bool)
        motif           = request.args.get('motif', default=None)
        rotation        = request.args.get('rotation', default=None, type=bool)
        resource = ShortTandemRepeat()
        if not all(v is None for v in [pathogenicity, motif, rotation]) and all(v is None for v in [uids, chrom, start, end, sources, location]):
            print(pathogenicity, file=sys.stdout)
            print(motif, file=sys.stdout)
            print(rotation, file=sys.stdout)
            result = short_tandem_repeat_queries(session, resource, pathogenicity, motif, rotation, limit, offset)
    else:
        raise BadRequest('valid resources are genes, short_tandem_repeats,\
             copy_number_variants, clinvar, conservation')
    # general queries 
    if result is None:
        result = genomic_feature_queries(session, resource, uids,chrom, start, end, sources, location, limit, offset)
    
    shutdown_session(session)
    # pass to create page
    result = create_page(resource, result, page, request.url)
    return jsonify(result)

# examples
# genes
# http://127.0.0.1:5000/api/v1/genesymbols
# http://127.0.0.1:5000/api/v1/genes?uids=1
# http://127.0.0.1:5000/api/v1/genes?names=LOC102725121
# http://127.0.0.1:5000/api/v1/genes?chrom=chr1&start=11868&end=14362
# http://127.0.0.1:5000/api/v1/genes?sources=refGene
# str 
# http://127.0.0.1:5000/api/v1/short_tandem_repeats?pathogenicity=True