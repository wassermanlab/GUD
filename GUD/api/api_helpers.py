from flask import request
from GUD import GUDUtils
from werkzeug.exceptions import NotFound, BadRequest
import math
import re
# print(names, file=sys.stdout)

## HELPER FUNCTIONS ##


def table_exists(table_name):
    engine = GUDUtils.get_engine()
    if not engine.dialect.has_table(engine, table_name): 
        raise BadRequest(table_name + ' table does not exist')

def get_db(db):
    if db == "hg19":
        GUDUtils.db = "hg19"
    elif db == "hg38":
        GUDUtils.db = "hg38"
    elif db == "test":
        GUDUtils.db = "test"
    else:
        raise BadRequest('database must be hg19 or hg38 or test')


def get_genomic_feature_results(resource, query, page_size, page) -> tuple:
    offset = (page-1)*page_size
    results = query.offset(offset).limit(page_size)
    results = [resource.as_genomic_feature(e) for e in results]
    results = [e.serialize() for e in results]
    return (query.count(), results)


def create_page(result_tuple, page, page_size, url) -> dict:
    """
    returns 404 error or a page
    """
    offset = (page-1)*page_size
    json = {}
    result_size = result_tuple[0]
    results = result_tuple[1]
    if result_size == 0:
        raise NotFound('No results from this query')
    if offset > result_size:
        raise BadRequest('Page range is invalid, valid range for this query is between 1 and ' +
                         str(math.ceil(result_size/page_size)))

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
    keys = get_mixin1_keys(request)

    # location query
    if (keys['start'] is not None and keys['end'] is not None and keys['location']
            is not None and keys['chrom'] is not None): 
        q = resource.select_by_location(
            session, keys['chrom'], keys['start'], keys['end'], keys['location'])
    else: 
        q = None
    # uid query
    if keys['uids'] is not None:
        q = resource.select_by_uids(session, q, keys['uids'])
    # sources query
    if keys['sources'] is not None:
        q = resource.select_by_sources(session, q, keys['sources'])
    return q


def genomic_feature_mixin2_queries(session, resource, request, query):
    keys = get_mixin2_keys(request)
    q = query
    if keys['experiments'] is not None:
        q = resource.select_by_experiments(session, q, keys['experiments'])
    if keys['samples'] is not None:
        q = resource.select_by_samples(session, q, keys['samples'])
    return q


def check_page(request):
    try:
        page = request.args.get('page', default=1, type=int)
    except:
        raise BadRequest('pages must be positive integers')
    if (page <= 0):
        raise BadRequest('pages must be positive integers')
    return page


def check_split(str_list, integer=False):
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

        if re.fullmatch('^chr(X|Y|[1-9]|1[0-9]|2[0-2])$', keys['chrom']) == None:
            raise BadRequest(
                "chromosome should be formatted as chrZ where z is X, Y, or 1-22")

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
