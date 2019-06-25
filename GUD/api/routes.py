from GUD.api import app
from GUD.api.db import establish_GUD_session, shutdown_session
from flask import request, jsonify
from GUD.ORM import Gene
from GUD.ORM.genomic_feature import GenomicFeature
import re
from werkzeug.exceptions import HTTPException, NotFound, BadRequest
import sys
# print(names, file=sys.stdout)
# errors
# 400 bad request, client input validation fails
# 404 not found


@app.route('/')
def index():
    return 'HOME'


def create_page(result: list, page, url) -> dict:
    """
    returns 404 error or a page
    """
    page_size = 20
    result_size = len(result)
    json = {}
    if page <= 0 or (page-1)*page_size > result_size:
        return 404  # invalid page
    start = (page-1)*page_size
    json = {'size': result_size,
            'results': result[start:start+page_size]}
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


@app.route('/api/v1/genesymbols')
def gene_symbols():
    url = request.url
    session = establish_GUD_session()
    page = int(request.args.get('page', default=1))
    result = Gene().get_all_gene_symbols(session)

    shutdown_session(session)
    result = create_page(result, page, url)
    if result == 404:
        raise NotFound('page range is invalid')
    return jsonify(result)


@app.route('/api/v1/<resource>')                       
def resource(resource):
    session = establish_GUD_session()
     # parameters
    page        = request.args.get('page', default=1, type=int)
    uids        = request.args.get('uids', default=None)
    chrom       = request.args.get('chrom', default=None)
    start       = request.args.get('start', default=None)
    end         = request.args.get('end', default=None)
    sources     = request.args.get('sources', default=None)
    location    = request.args.get('location', default='within', type=str)
    if (resource == 'genes'):
        names       = request.args.get('names', default=None)
        resource    = Gene()
        if names is not None:
            names = names.split(',')
            result = resource.select_by_names(session, names)
    elif (resource == 'short_tandem_repeats'):
        pass
    else: 
        raise BadRequest('valid resources are genes, short_tandem_repeats,\
             copy_number_variants, clinvar, conservation')  
    
    if uids is not None:
        try:
            uids = uids.split(',')
            uids = [int(e) for e in uids]
        except:
            raise BadRequest("uids must be positive integers seperated by commas (,).")
        result = resource.select_by_uids(session, uids)
    elif chrom is not None and start is not None and end is not None:
        try: 
            start = int(start) - 1
            end = int(end)
        except:
            raise BadRequest("start and end should be formatted as integers, \
            chromosomes should be formatted as chrZ.")
        if location == 'exact':
            result = resource.select_by_exact_location(session, chrom, start, end)
        else:
            result = resource.select_by_location(session, chrom, start, end)

    shutdown_session(session)
    if len(result) == 0:                                                        # 404 if nothing is found
        raise NotFound("no refgene genes found by the search.")     
    result = [resource.as_genomic_feature(e) for e in result]                   # turn to genomicFeature
    result = [e.serialize() for e in result]                                    # serialize to json
    result = create_page(result, page, request.url)                                     # pass to create page
    if result == 404:                                                           # if page range is incorrect    
        raise NotFound('page range is invalid')     
    return jsonify(result)
        

# @app.route('/api/v1/genes')
# def genes():
#     session = establish_GUD_session()
#     # parameters
#     page        = request.args.get('page', default=1, type=int)
#     names       = request.args.get('names', default=None)
#     uids        = request.args.get('uids', default=None)
#     chrom       = request.args.get('chrom', default=None)
#     start       = request.args.get('start', default=None)
#     end         = request.args.get('end', default=None)
#     sources     = request.args.get('sources', default=None)
#     location    = request.args.get('location', default="exact") ## options, within or exact 
#     # queries
#     ## TODO: add pagination 
#     gene = Gene()
#     if uids is not None:
#         try:
#             uids = uids.split(',')
#             uids = [int(e) for e in uids]
#         except:
#             raise BadRequest("uids must be positive integers seperated by commas (,).")
#         result = gene.select_by_uids(session, uids)
#     elif names is not None:
#         names = names.split(',')
#         result = gene.select_by_names(session, names)
#     else:
#         try: 
#             start = int(start) - 1
#             end = int(end)
#         except:
#             raise BadRequest("start and end should be formatted as integers, \
#             chromosomes should be formatted as chrZ.")
#         if location == 'exact':
#             result = gene.select_by_exact_location(session, chrom, start, end)
#         else:
#             result = gene.select_by_location(session, chrom, start, end)
#     shutdown_session(session)
#     if len(result) == 0:
#         raise NotFound("no refgene genes found by the search.")
#     result = [gene.as_genomic_feature(e) for e in result]
#     result = [e.serialize() for e in result]
#     return jsonify(result)



# examples
# http://127.0.0.1:5000/api/v1/genesymbols
# http://127.0.0.1:5000/api/v1/genes?uids=1
# http://127.0.0.1:5000/api/v1/genes?names=LOC102725121
# http://127.0.0.1:5000/api/v1/genes?chrom=chr1&start=11868&end=14362