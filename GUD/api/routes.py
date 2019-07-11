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

## HELPER FUNCTIONS ##


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

    if resource is not False:
        results = [resource.as_genomic_feature(e) for e in results]
        results = [e.serialize() for e in results]
    else:
        results = [g[0] for g in results]

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

    q = resource.select_by_location(
        session, keys['chrom'], keys['start'], keys['end'], keys['location'])
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
    clinvarIDs = check_split(request.args.get(
        'clinvar_ids', default=None), True)
    q = genomic_feature_mixin1_queries(session, clinvar, request)
    if clinvarIDs is not None:
        q = clinvar.select_by_clinvarID(session, q, clinvarIDs)
    shutdown_session(session)
    result = create_page(clinvar, q, page, request.url)
    return jsonify(result)


@app.route('/api/v1/copy_number_variants')
def mixin1():
    page = check_page(request)
    session = establish_GUD_session()
    resource = CNV()
    q = genomic_feature_mixin1_queries(session, resource, request)
    shutdown_session(session)
    result = create_page(resource, q, page, request.url)
    return jsonify(result)


@app.route('/api/v1/genes')
def genes():
    page = check_page(request)
    session = establish_GUD_session()
    resource = Gene()
    q = genomic_feature_mixin1_queries(session, resource, request)
    names = check_split(request.args.get('names', default=None))
    if names is not None:
        q = resource.select_by_names(session, q, names)
    shutdown_session(session)
    result = create_page(resource, q, page, request.url)
    return jsonify(result)


@app.route('/api/v1/genes/symbols')
def gene_symbols():
    url = request.url
    session = establish_GUD_session()
    page = int(request.args.get('page', default=1))
    q = Gene().get_all_gene_symbols(session)
    shutdown_session(session)
    result = create_page(False, q, page, url)
    return jsonify(result)


@app.route('/api/v1/short_tandem_repeats')
def strs():
    page = check_page(request)
    session = establish_GUD_session()
    resource = ShortTandemRepeat()
    q = genomic_feature_mixin1_queries(session, resource, request)
    names = check_split(request.args.get('names', default=None))
    try:
        motif = request.args.get('motif', default=None, type=str)
        rotation = request.args.get('rotation', default=False, type=bool)
    except:
        raise BadRequest(
            'rotation must be set to True or False if motif is given')
    if motif is not None:
        q = resource.select_by_motif(session, motif, q, rotation)
    result = create_page(resource, q, page, request.url)
    return jsonify(result)


@app.route('/api/v1/short_tandem_repeats/pathogenic')
def pathogenic_strs():
    page = check_page(request)
    session = establish_GUD_session()
    resource = ShortTandemRepeat()
    q = resource.select_by_pathogenicity(session)
    result = create_page(resource, q, page, request.url)
    return jsonify(result)


@app.route('/api/v1/enhancers')
@app.route('/api/v1/dna_accessibility')
def mixin2():
    page = check_page(request)
    session = establish_GUD_session()
    if (request.path == '/api/v1/dna_accessibility'):
        resource = DNAAccessibility()
    elif (request.path == '/api/v1/enhancers'):
        resource = Enhancer()
    q = genomic_feature_mixin1_queries(session, resource, request)
    q = genomic_feature_mixin2_queries(session, resource, request, q)
    shutdown_session(session)
    result = create_page(resource, q, page, request.url)
    return jsonify(result)

@app.route('/api/v1/histone_modifications')
def histone_modifications():
    page = check_page(request)
    session = establish_GUD_session()
    resource = HistoneModification()
    histone_types = check_split(request.args.get('histone_types', default=None))
    q = genomic_feature_mixin1_queries(session, resource, request)
    q = genomic_feature_mixin2_queries(session, resource, request, q)
    if histone_types is not None: 
        q = resource.select_by_histone_type(q, histone_types)
    shutdown_session(session) 
    result = create_page(resource, q, page, request.url)
    return jsonify(result)

@app.route('/api/v1/tads')
def tads():
    page = check_page(request)
    session = establish_GUD_session()
    resource = TAD()
    restriction_enzymes = check_split(request.args.get('restriction_enzymes', default=None))
    q = genomic_feature_mixin1_queries(session, resource, request)
    q = genomic_feature_mixin2_queries(session, resource, request, q)
    if restriction_enzymes is not None:
        q = resource.select_by_restriction_enzymes(q, restriction_enzymes)
    shutdown_session(session) 
    result = create_page(resource, q, page, request.url)
    return jsonify(result)

@app.route('/api/v1/tf_binding')
def tf_binding():
    page = check_page(request)
    session = establish_GUD_session()
    resource = TFBinding()
    tfs = check_split(request.args.get('tfs', default=None))
    q = genomic_feature_mixin1_queries(session, resource, request)
    q = genomic_feature_mixin2_queries(session, resource, request, q)
    if tfs is not None:
        q = resource.select_by_tf(q, tfs)
    shutdown_session(session) 
    result = create_page(resource, q, page, request.url)
    return jsonify(result)

@app.route('/api/v1/tss')
def tss():
    samples = request.args.get('samples', default=None, type=str)
    if samples is not None:
        return BadRequest('Cannot query TSS table by sample')
    page = check_page(request)
    session = establish_GUD_session()
    resource = TSS()
    genes = check_split(request.args.get('genes', default=None))
    q = genomic_feature_mixin1_queries(session, resource, request)
    q = genomic_feature_mixin2_queries(session, resource, request, q)
    if genes is not None: 
        q = resource.select_by_genes(q, genes)
    shutdown_session(session) 
    result = create_page(resource, q, page, request.url)
    return jsonify(result)

@app.route('/api/v1/tss/genic')
def genic_tss():
    page = check_page(request)
    session = establish_GUD_session()
    resource = TSS()
    q = resource.select_all_genic_tss(session)
    result = create_page(resource, q, page, request.url)
    return jsonify(result)

@app.route('/')
def index():
    return 'HOME'
# print(names, file=sys.stdout)