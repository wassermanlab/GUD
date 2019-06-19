from GUD.api import app
from GUD.api.db import establish_GUD_session, shutdown_session 
from flask import request, jsonify
from GUD.ORM import Gene
from GUD.ORM.genomic_feature import GenomicFeature
import sys

@app.route('/')
def index():
    return 'HOME'

@app.route('/genes')
def genes():
    session = establish_GUD_session()

    # parameters
    regionID    = request.args.get('regionID', default=None)
    sourceID    = request.args.get('sourceID', default=None)
    names       = request.args.get('names', default=None)
    uid         = request.args.get('uid', default=None)
    chrom       = request.args.get('chrom', default=None)
    start       = request.args.get('start', default=None)
    end         = request.args.get('end', default=None)
    # queries 
    print(chrom, file=sys.stdout)
    print(start, file=sys.stdout)
    print(end, file=sys.stdout)
    gene = Gene()
    if uid is not None:
       result = gene.select_by_uid(session, uid, True) 
    elif names is not None: 
        names = names.split(',')
        # print(names, file=sys.stdout)
        result = gene.select_by_names(session, names, True)
    else:
        result = gene.select_by_location(session, chrom, int(start), int(end), True)
    shutdown_session(session)
    if type(result) is list:
        result = [e.serialize() for e in result]
    else: 
        result = result.serialize()
    return jsonify(result)

## examples
#http://127.0.0.1:5000/genes?uid=1
#http://127.0.0.1:5000/genes?names=LOC102725121
#http://127.0.0.1:5000/genes?chrom=chr1&start=11868&end=14362