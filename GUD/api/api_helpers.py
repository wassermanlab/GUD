import time
from flask import request, jsonify
from GUD import GUDUtils
from werkzeug.exceptions import NotFound, BadRequest
import math
import re
from sqlalchemy import func
from GUD.ORM import ShortTandemRepeat
import time

## HELPER FUNCTIONS ##
def get_result_from_query(query, request, resource, result_tuple_type="simple"):
    if query is None:
        raise BadRequest('query not specified correctly')
    # results = query.all() 
    results = []
    start = time.time()
    for q in query:
        s = time.time()
        r = q.all()
        print(q.statement.compile(compile_kwargs={"literal_binds": True}))
        results = results + r
        print('getting subquery took', time.time()-s, 'seconds.')
    print('getting queries took', time.time()-start, 'seconds.')
    # print(results.statement.compile(compile_kwargs={"literal_binds": True}))
    # serialize and get uids of first and last element returned
    start = time.time()
    if (result_tuple_type == "genomic_feature"):
        results = [resource.as_genomic_feature(e) for e in results]
    results = [e.serialize() for e in results]
    results = create_page(results, request.url)
    print('serializing', time.time()-start, 'seconds.')
    return jsonify(results)


def create_page(results, url) -> dict:
    """
    returns 404 error or a page
    """
    json = {}
    if len(results) == 0:
        raise NotFound('No results from this query')
    json = {'results': results}
    return json


def table_exists(table_name, engine):
    if not engine.dialect.has_table(engine, table_name):
        raise BadRequest(table_name + ' table does not exist')


def set_db(db):
    if db == "grch37":
        GUDUtils.db = "grch37"
    elif db == "grch38":
        GUDUtils.db = "grch38"
    else:
        raise BadRequest(
            'database must be grch37 or grch38')


def genomic_feature_mixin1_queries(session, resource, request):
    """make genomic feature 1 queries"""
    keys = get_mixin1_keys(request)
    # query array 
    q = []
    # all location
    if (keys['merged_start_end'] is not None and keys['location'] is not None and keys['chrom'] is not None):
        for i in keys['merged_start_end']:
            q.append(resource.select_by_location(session, None, keys['chrom'], i[0], i[1], keys['location']))
    else: 
        q.append(resource.select_all(session, None))
    if keys['uids'] is not None:
        q = [resource.select_by_uids(session, i, keys['uids']) for i in q]
    # sources query
    if keys['sources'] is not None:
        q = [resource.select_by_sources(session, i, keys['sources']) for i in q]
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
            'sources': [], 
            'last_uid': '',
            'merged_start_end': None}
            
    keys['chrom'] = request.args.get('chrom', default=None, type=str)
    keys['end'] = check_split(request.args.get('end', default=None))
    keys['location'] = request.args.get('location', default=None, type=str)
    keys['start'] = check_split(request.args.get('start', default=None))
    keys['sources'] = check_split(request.args.get('sources', default=None))
    keys['uids'] = check_split(request.args.get('uids', default=None))
    keys['last_uid'] = request.args.get('last_uid', default=0, type=int)

    if keys['uids'] is not None:        # convert uids if they are in uri
        for i in range(len(keys['uids'])):
            if keys['uids'][i].isdigit():
                keys['uids'][i] = int(keys['uids'][i])
    
    # check that location is specified
    

    if (keys['start'] is not None and keys['end'] is not None and keys['location']
            is not None and keys['chrom'] is not None):
        if (len(keys['start']) != len(keys['end'])):
            raise BadRequest("start and end lists should be the same length")
        try:
            keys['merged_start_end'] = [(int(keys['start'][i]), int(keys['end'][i])) for i in range(0, len(keys['start']))] 
        except:
            raise BadRequest("start and end should be formatted as integers")
        for i in keys['merged_start_end']:
            if ((i[1]-i[0]) > 4000000): # check limit 
                raise BadRequest("each region must be less than 4,000,000bp")
        if re.fullmatch('^(X|Y|[1-9]|1[0-9]|2[0-2])$', keys['chrom']) == None:
            raise BadRequest(
                "chromosome should be formatted as Z where Z is X, Y, or 1-22")
        if keys['location'] not in ['within', 'overlapping', 'exact']:
            raise BadRequest(
                "location must be specified as within, overlapping, or exact")
    return keys
