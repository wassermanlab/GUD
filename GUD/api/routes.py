from GUD.api import app
from GUD.api.db import establish_GUD_session, shutdown_session
from flask import request, jsonify
from GUD.ORM import (Gene, ShortTandemRepeat, CNV, ClinVar, Conservation,
                     DNAAccessibility, Enhancer, HistoneModification, TAD, TFBinding, TSS)
from GUD.ORM.genomic_feature import GenomicFeature
from werkzeug.exceptions import HTTPException, NotFound, BadRequest
import sys, math, re
# print(names, file=sys.stdout)


def create_page(resource, query, page, url) -> dict:
    """
    returns 404 error or a page
    """
    page_size = 20
    offset = (page-1)*page_size
    json = {}
    result_size = query.count()
    results = query.offset(offset).limit(page_size)
    if result_size == 0:
        raise NotFound('No results from this query')
    if offset > result_size:
        raise BadRequest('Page range is invalid, valid range for this query is between 1 and ' +
                         str(math.ceil(result_size/page_size)-1))

    if (resource != False):
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


def genomic_feature_mixin1_queries(session, resource, request):
    keys = set(request.args)
    keys.discard('page')
    keys = get_mixin1_keys(request)

    # location query
    q = resource.select_by_location(
        session, keys['chrom'], keys['start'], keys['end'], keys['location'])
    # add uids
    if keys['uids'] is not None:
        # add sources
        q = resource.select_by_uids(session, q, keys['uids'])
    if keys['sources'] is not None:
        q = resource.select_by_sources(session, q, keys['sources'])
    return q

## HELPER FUNCTIONS ##


def check_page(request):
    try:
        page = request.args.get('page', default=1, type=int)
    except:
        raise BadRequest('pages must be positive integers')
    if (page <= 0):
        raise BadRequest('pages must be positive integers')
    return page


def check_split(str_list, integer = False):
    if str_list == None:
        return None
    s = str_list.split(',')
    if len(str_list) > 1000 or len(str_list) < 1:
        raise BadRequest(
            "list query parameters must be greater than 0 and less than 1000")
    if integer is True:
        try:
            s = [int(i) for i in s]
        except:
            raise BadRequest('list query parameter needs to be integer')
    
    return s


def get_mixin1_keys(request):
    keys = {'chrom': '',
            'start': '',
            'end': '',
            'location': '',
            'sources': []}
    keys['chrom'] = request.args.get('chrom', default=None, type=str)
    keys['end'] = request.args.get('end', default=None)
    keys['location'] = request.args.get('location', default=None, type=str)
    keys['start'] = request.args.get('start', default=None)
    keys['sources'] = check_split(request.args.get('sources', default=None))
    keys['uids'] = check_split(request.args.get('uids', default=None))

    if keys['uids'] is not None:        # convert uids if they are in uri
        for i in range(len(keys['uids'])):
            if keys['uids'][i].isdigit():
                keys['uids'][i] = int(keys['uids'][i])

    try:
        keys['start'] = int(keys['start']) - 1
        keys['end'] = int(keys['end'])
    except:
        raise BadRequest("start and end should be formatted as integers")

    if re.fullmatch('^chr(X|Y|[1-9]|1[0-9]|2[0-2])$', keys['chrom']) == None:
        raise BadRequest("chromosome should be formatted as chrZ where z \
            is X, Y, or 1-22")

    if (keys['end'] - keys['start']) > 4000000:
        raise BadRequest("start and end must be less than 4000000bp apart")

    if keys['location'] not in ['within', 'overlapping', 'exact']:
        raise BadRequest(
            "location must be specified as withing, overlapping, or exact")
    return keys


def get_mixin2_keys(request):
    keys = {'experiments': [],
            'samples': []}
    keys['experiments'] = check_split(
        request.args.get('experiments', default=None))
    keys['samples'] = check_split(request.args.get('samples', default=None))
    return keys


## ROUTES ##

@app.route('/api/v1/clinvar')
def clinvar():
    page = check_page(request)
    session = establish_GUD_session()
    clinvar = ClinVar()
    clinvarIDs = check_split(request.args.get('clinvar_ids', default=None), True)
    q = genomic_feature_mixin1_queries(session, clinvar, request)
    if clinvarIDs is not None:
        q = clinvar.select_by_clinvarID(session, q, clinvarIDs)
    shutdown_session(session)
    result = create_page(clinvar, q, page, request.url)
    return jsonify(result)


@app.route('/')
def index():
    return 'HOME'
