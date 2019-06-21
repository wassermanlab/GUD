from GUD.api import app
from GUD.api.db import establish_GUD_session, shutdown_session
from flask import request, jsonify
from GUD.ORM import Gene
from GUD.ORM.genomic_feature import GenomicFeature
import re
from werkzeug.exceptions import HTTPException, NotFound
import sys
# print(names, file=sys.stdout)


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


# @app.route('/api/v1/genes')
# def genes():
#     page_size = 20
#     session = establish_GUD_session()
#     # parameters
#     page        = request.args.get('page', default=None)
#     names       = request.args.get('names', default=None)
#     uids         = request.args.get('uids', default=None)
#     chrom       = request.args.get('chrom', default=None)
#     start       = request.args.get('start', default=None)
#     end         = request.args.get('end', default=None)
#     sources     = request.args.get('sources', default=None)
#     # queries
#     gene = Gene()
#     if uids is not None:
#         uids = uids.split(',')
#         uids = map(int, uids)
#         result = gene.select_by_uids(session, uids)
#     elif names is not None:
#         names = names.split(',')
#         # print(names, file=sys.stdout)
#         result = gene.select_by_names(session, names, True)
#     else:
#         result = gene.select_by_location(session, chrom, int(start), int(end), True)
#     shutdown_session(session)
#     if type(result) is list:
#         result = [e.serialize() for e in result]
#     else:
#         result = result.serialize()
#     return jsonify(result)


# examples
# http://127.0.0.1:5000/api/v1/genesymbols
# http://127.0.0.1:5000/api/v1/genes?uid=1
# http://127.0.0.1:5000/api/v1/genes?names=LOC102725121
# http://127.0.0.1:5000/api/v1/genes?chrom=chr1&start=11868&end=14362
