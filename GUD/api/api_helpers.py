from flask import request, jsonify
from GUD import GUDUtils
from werkzeug.exceptions import NotFound, BadRequest
import math
import re
from sqlalchemy import func
from GUD.ORM import ShortTandemRepeat
import time

## HELPER FUNCTIONS ##
def get_result_from_query(query, request, resource, page_size=20, result_tuple_type="simple"):
    last_uid = request.args.get('last_uid', default=0, type=int)
    if query is None:
        raise BadRequest('query not specified correctly')
    results = query.filter(type(resource).uid > last_uid).with_hint(type(resource), 'USE INDEX (PRIMARY)')\
        .order_by(type(resource).uid).limit(page_size)  # seek method for paginating, we must specify the index to have speed up
    print(results)   # TODO: take this out later
    # serialize and get uids of first and last element returned
    try:
        if (result_tuple_type == "genomic_feature"):
            last_uid = getattr(results[page_size-1], type(resource).__name__).uid
        else:
            last_uid = results[page_size-1].uid
    except:
        last_uid = None
    if (result_tuple_type == "genomic_feature"):
        results = [resource.as_genomic_feature(e) for e in results]
    results = [e.serialize() for e in results]
    results = create_page(results, last_uid, page_size, request.url)
    return jsonify(results)


def create_page(results, last_uid, page_size, url) -> dict:
    """
    returns 404 error or a page
    """
    json = {}
    if len(results) == 0:
        raise NotFound('No results from this query')
    json = {'results': results}
    if last_uid != None: 
        if (re.search('\?', url) is None):
            next_page = url+'?last_uid='+str(last_uid)
        elif (re.search('last_uid', url) is None):
            next_page = url+'&last_uid='+str(last_uid)
        else:
            next_page = re.sub('last_uid=\d+', 'last_uid='+str(last_uid), url)
        json['next'] = next_page
    return json


def table_exists(table_name, engine):
    if not engine.dialect.has_table(engine, table_name):
        raise BadRequest(table_name + ' table does not exist')


def set_db(db):
    if db == "hg19":
        GUDUtils.db = "hg19"
    elif db == "hg38":
        GUDUtils.db = "hg38"
    elif db == "test":
        GUDUtils.db = "test"
    elif db == "test_hg38_chr22":
        GUDUtils.db = "test_hg38_chr22"
    else:
        raise BadRequest(
            'database must be hg19 or hg38 or test or test_hg38_chr22')


def genomic_feature_mixin1_queries(session, resource, request):
    """make genomic feature 1 queries"""
    keys = get_mixin1_keys(request)
    # location query
    q = resource.select_all(session, None)
    # just chrom
    if (keys['start'] is None and keys['end'] is None and keys['location'] is None and keys['chrom'] is not None):
        q = resource.select_by_location(session, q, keys['chrom'])
    # all location
    elif (keys['start'] is not None and keys['end'] is not None and keys['location'] is not None and keys['chrom'] is not None):
        q = resource.select_by_location(
                session, q, keys['chrom'], keys['start'], keys['end'], keys['location'])
    # partial location 
    elif (keys['start'] is not None or keys['end'] is not None or keys['location'] is not None or keys["chrom"] is not None):
        raise BadRequest("To filter by location you must specify location, chrom, start, and end or just a chrom.")
    # uid query
    if keys['uids'] is not None:
        q = resource.select_by_uids(session, q, keys['uids'])
    # sources query
    if keys['sources'] is not None:
        q = resource.select_by_sources(session, q, keys['sources'])
    return q


def genomic_feature_mixin2_queries(session, resource, request, query):
    """make genomic feature 2 queries"""
    keys = get_mixin2_keys(request)
    q = query
    if keys['experiments'] is not None:
        q = resource.select_by_experiments(session, q, keys['experiments'])
    if keys['samples'] is not None:
        q = resource.select_by_samples(session, q, keys['samples'])
    return q


def check_split(str_list, integer=False):
    """split string delimeted by ','"""
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

    if (keys['start'] is not None and keys['end'] is not None and keys['location']
            is not None and keys['chrom'] is not None):
        try:
            keys['start'] = int(keys['start']) - 1
            keys['end'] = int(keys['end'])
        except:
            raise BadRequest("start and end should be formatted as integers")
        if re.fullmatch('^(X|Y|[1-9]|1[0-9]|2[0-2])$', keys['chrom']) == None:
            raise BadRequest(
                "chromosome should be formatted as Z where Z is X, Y, or 1-22")
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
