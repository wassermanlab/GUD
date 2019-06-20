from GUD.api import app
from GUD.api.db import establish_GUD_session, shutdown_session 
from flask import request, jsonify
from GUD.ORM import Gene
from GUD.ORM.genomic_feature import GenomicFeature
import sys


@app.route('/')
def index():
    return 'HOME'

@app.route('/api/v1/genesymbols')
def gene_symbols():
    page_size   = 40
    session     = establish_GUD_session()
    page        = int(request.args.get('page', default=1))
    result = Gene().get_all_gene_symbols(session)
    shutdown_session(session)
    result_size = len(result)

    if page == 1: 
    ## if page is 1 
    ## if page is within range
    ## if page is not within range -> 404 error i t hink
        start = (page-1)*40 
        json = {
                'next': '', 
                'size': result_size,
                'results': result[start:40]}
    return jsonify(json)


@app.route('/api/v1/genes')
def genes():
    page_size = 20
    session = establish_GUD_session()

    # parameters
    page        = request.args.get('page', default=None)
    regionID    = request.args.get('regionID', default=None)
    sourceID    = request.args.get('sourceID', default=None)
    names       = request.args.get('names', default=None)
    uids         = request.args.get('uids', default=None)
    chrom       = request.args.get('chrom', default=None)
    start       = request.args.get('start', default=None)
    end         = request.args.get('end', default=None)
    # queries 
    gene = Gene()
    if uids is not None:
        uids = uids.split(',')
        uids = map(int, uids)
        result = gene.select_by_uids(session, uids, True) 
    # elif names is not None: 
    #     names = names.split(',')
    #     # print(names, file=sys.stdout)
    #     result = gene.select_by_names(session, names, True)
    # else:
    #     result = gene.select_by_location(session, chrom, int(start), int(end), True)
    shutdown_session(session)
    if type(result) is list:
        result = [e.serialize() for e in result]
    else: 
        result = result.serialize()
    return jsonify(result)

## examples
#http://127.0.0.1:5000/api/v1/genesymbols
#http://127.0.0.1:5000/api/v1/genes?uid=1
#http://127.0.0.1:5000/api/v1/genes?names=LOC102725121
#http://127.0.0.1:5000/api/v1/genes?chrom=chr1&start=11868&end=14362