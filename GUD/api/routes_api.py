from GUD import GUDUtils
from GUD.api import app
from flask import request, jsonify
from GUD.ORM import (Gene, ShortTandemRepeat, CNV, ClinVar, Conservation,
                     DNAAccessibility, Enhancer, HistoneModification, TAD, 
                     TFBinding, TSS, Chrom, Sample, Experiment, Source, Expression)
from GUD.api.api_helpers import *
from werkzeug.exceptions import HTTPException, NotFound, BadRequest
import sys
page_size = 20
# print(names, file=sys.stdout)

@app.route('/api/v1/<db>/clinvar')
def clinvar(db):

    get_db(db)
    db_name = GUDUtils._get_db_name()
    engine, Session = GUDUtils._get_engine_session(db_name)

    table_exists('clinvar')
    page = check_page(request)
    resource = ClinVar()
    clinvarIDs = check_split(request.args.get(
        'clinvar_ids', default=None), True)
    q = genomic_feature_mixin1_queries(Session, resource, request)
    if clinvarIDs is not None:
        q = resource.select_by_clinvarID(Session, q, clinvarIDs)
    
    Session.close()
    engine.dispose()

    result_tuple = get_genomic_feature_results(resource, q, page_size,  page)
    result = create_page(result_tuple, page, page_size, request.url)
    response = jsonify(result)
    return response

@app.route('/api/v1/<db>/copy_number_variants')
def copy_number_variants(db):
    get_db(db)
    db_name = GUDUtils._get_db_name()
    engine, Session = GUDUtils._get_engine_session(db_name)

    table_exists('copy_number_variants')
    page = check_page(request)
    clinical_assertion  = request.args.get('clinical_assertion', default=None)
    clinvar_accession   = request.args.get('clinvar_accession', default=None)
    dbVar_accession     = request.args.get('dbVar_accession', default=None)
    resource = CNV()
    q = genomic_feature_mixin1_queries(Session, resource, request)
    
    if clinical_assertion is not None:
        q = resource.select_by_clinical_assertion(Session, q, clinical_assertion)
    if clinvar_accession is not None:
        q = resource.select_by_clinvar_accession(Session, q, clinvar_accession)
    if dbVar_accession is not None:
        q = resource.select_by_dbvar_accession(Session, q, dbVar_accession)
    
    Session.close()
    engine.dispose()
    
    result_tuple = get_genomic_feature_results(resource, q, page_size,  page)
    result = create_page(result_tuple, page, page_size, request.url)
    response = jsonify(result)
    return response

@app.route('/api/v1/<db>/genes')
def genes(db):
    get_db(db)
    db_name = GUDUtils._get_db_name()
    engine, Session = GUDUtils._get_engine_session(db_name)

    table_exists('genes')
    page = check_page(request)
    
    resource = Gene()
    q = genomic_feature_mixin1_queries(Session, resource, request)
    names = check_split(request.args.get('names', default=None))
    if names is not None:
        q = resource.select_by_names(Session, q, names)
    
    Session.close()
    engine.dispose()

    if (q is None):
        raise BadRequest('Query not specified correctly')
    result_tuple = get_genomic_feature_results(resource, q, page_size,  page)
    result = create_page(result_tuple, page,page_size, request.url)
    response = jsonify(result)
    return response

@app.route('/api/v1/<db>/genes/symbols')
def gene_symbols(db):
    get_db(db)
    db_name = GUDUtils._get_db_name()
    engine, Session = GUDUtils._get_engine_session(db_name)
    table_exists('genes')
    
    url = request.url
    
    page = int(request.args.get('page', default=1))
    q = Gene().get_all_gene_symbols(Session)
    
    Session.close()
    engine.dispose()
    page_size = 1000 ## set custom page size 
    offset = (page-1)*page_size
    result = q.offset(offset).limit(page_size)
    result =  [g[0] for g in result]
    result_tuple = (q.count(), result)
    result = create_page(result_tuple, page, page_size, request.url)                                 
    page_size = 20  ## reset page size
    response = jsonify(result)
    return response

@app.route('/api/v1/<db>/short_tandem_repeats')
def strs(db):
    get_db(db)
    db_name = GUDUtils._get_db_name()
    engine, Session = GUDUtils._get_engine_session(db_name)
    table_exists('short_tandem_repeats')
    page = check_page(request)
    
    resource = ShortTandemRepeat()
    q = genomic_feature_mixin1_queries(Session, resource, request)
    try:
        motif = request.args.get('motif', default=None, type=str)
        rotation = request.args.get('rotation', default=False, type=bool)
    except:
        raise BadRequest(
            'rotation must be set to True or False if motif is given')
    if motif is not None:
        q = resource.select_by_motif(Session, motif, q, rotation)
    result_tuple = get_genomic_feature_results(resource, q, page_size,  page)
    result = create_page(result_tuple, page, page_size, request.url) 
    Session.close()
    engine.dispose()
    response = jsonify(result)
    return response

@app.route('/api/v1/<db>/short_tandem_repeats/pathogenic')
def pathogenic_strs(db):
    get_db(db)
    db_name = GUDUtils._get_db_name()
    engine, Session = GUDUtils._get_engine_session(db_name)

    table_exists('short_tandem_repeats')
    page = check_page(request)
    
    resource = ShortTandemRepeat()
    q = resource.select_by_pathogenicity(Session)
    result_tuple = get_genomic_feature_results(resource, q, page_size,  page)
    result = create_page(result_tuple, page, page_size, request.url)
    
    Session.close()
    engine.dispose()
    return jsonify(result)

@app.route('/api/v1/<db>/enhancers')
@app.route('/api/v1/<db>/dna_accessibility')
def mixin2(db):
    get_db(db)
    db_name = GUDUtils._get_db_name()
    engine, Session = GUDUtils._get_engine_session(db_name)
  
    page = check_page(request)
    
    if ('dna_accessibility' in request.path):
        table_exists('dna_accessibility')
        resource = DNAAccessibility()
    elif ('enhancers' in request.path):
        table_exists('enhancers')
        resource = Enhancer()
    q = genomic_feature_mixin1_queries(Session, resource, request)
    q = genomic_feature_mixin2_queries(Session, resource, request, q)
    
    Session.close()
    engine.dispose()

    result_tuple = get_genomic_feature_results(resource, q, page_size,  page)
    result = create_page(result_tuple, page, page_size, request.url)
    return jsonify(result)

@app.route('/api/v1/<db>/histone_modifications')
def histone_modifications(db):
    get_db(db)
    db_name = GUDUtils._get_db_name()
    engine, Session = GUDUtils._get_engine_session(db_name)
 
    table_exists('histone_modifications')
    page = check_page(request)

    resource = HistoneModification()
    histone_types = check_split(request.args.get('histone_types', default=None))
    q = genomic_feature_mixin1_queries(Session, resource, request)
    q = genomic_feature_mixin2_queries(Session, resource, request, q)

    if histone_types is not None: 
        q = resource.select_by_histone_type(q, histone_types)
    Session.close()
    engine.dispose()

    result_tuple = get_genomic_feature_results(resource, q, page_size,  page)
    result = create_page(result_tuple, page,page_size, request.url)
    response = jsonify(result)
    return response

@app.route('/api/v1/<db>/tads')
def tads(db):
    get_db(db)
    db_name = GUDUtils._get_db_name()
    engine, Session = GUDUtils._get_engine_session(db_name)

    engine.dispose()
    table_exists('tads')
    page = check_page(request)
    
    resource = TAD()
    restriction_enzymes = check_split(request.args.get('restriction_enzymes', default=None))
    q = genomic_feature_mixin1_queries(Session, resource, request)
    q = genomic_feature_mixin2_queries(Session, resource, request, q)
    if restriction_enzymes is not None:
        q = resource.select_by_restriction_enzymes(q, restriction_enzymes)
    
    Session.close()
    engine.dispose()

    result_tuple = get_genomic_feature_results(resource, q, page_size,  page)
    result = create_page(result_tuple, page, page_size, request.url)
    return jsonify(result)

@app.route('/api/v1/<db>/tf_binding')
def tf_binding(db):
    get_db(db)
    db_name = GUDUtils._get_db_name()
    engine, Session = GUDUtils._get_engine_session(db_name)

    table_exists('tf_binding')
    page = check_page(request)
    
    resource = TFBinding()
    tfs = check_split(request.args.get('tfs', default=None))
    q = genomic_feature_mixin1_queries(Session, resource, request)
    q = genomic_feature_mixin2_queries(Session, resource, request, q)
    if tfs is not None:
        q = resource.select_by_tf(q, tfs)
    
    Session.close()
    engine.dispose()
    result_tuple = get_genomic_feature_results(resource, q, page_size,  page)
    result = create_page(result_tuple, page, page_size, request.url)
    return jsonify(result)

@app.route('/api/v1/<db>/tss')
def tss(db):
    get_db(db)
    db_name = GUDUtils._get_db_name()
    engine, Session = GUDUtils._get_engine_session(db_name)

    table_exists('transcription_start_sites')
    samples = request.args.get('samples', default=None, type=str)
    if samples is not None:
        return BadRequest('Cannot query TSS table by sample')
    page = check_page(request)
    
    resource = TSS()
    genes = check_split(request.args.get('genes', default=None))
    q = genomic_feature_mixin1_queries(Session, resource, request)
    q = genomic_feature_mixin2_queries(Session, resource, request, q)
    if genes is not None: 
        q = resource.select_by_genes(q, genes)
    
    Session.close()
    engine.dispose()

    result_tuple = get_genomic_feature_results(resource, q, page_size,  page)
    result = create_page(result_tuple, page, page_size, request.url)
    return jsonify(result)

@app.route('/api/v1/<db>/tss/genic')
def genic_tss(db):
    get_db(db)
    db_name = GUDUtils._get_db_name()
    engine, Session = GUDUtils._get_engine_session(db_name)
    
    table_exists('transcription_start_sites')
    page = check_page(request)
    
    resource = TSS()
    q = resource.select_all_genic_tss(Session)
    Session.close()
    engine.dispose()
    result_tuple = get_genomic_feature_results(resource, q, page_size,  page)
    result = create_page(result_tuple, page, page_size, request.url)
    return jsonify(result)

@app.route('/api/v1/<db>/chroms')
def chroms(db):
    get_db(db)
    db_name = GUDUtils._get_db_name()
    engine, Session = GUDUtils._get_engine_session(db_name)
 
    table_exists('chroms')
    page = check_page(request)
    
    resource = Chrom()
    q = resource.select_all_chroms(Session)
    offset = (page-1)*page_size
    results = q.offset(offset).limit(page_size)
    results = [e.serialize() for e in results]
    result_tuple = (q.count(), results)
    result = create_page(result_tuple, page, page_size, request.url)
    
    Session.close()
    engine.dispose()
    return jsonify(result)

@app.route('/api/v1/<db>/sources')
def sources(db):
    get_db(db)
    db_name = GUDUtils._get_db_name()
    engine, Session = GUDUtils._get_engine_session(db_name)

    table_exists('sources')
    page = check_page(request)
    
    resource = Source()
    q = resource.select_all_sources(Session)
    offset = (page-1)*page_size
    results = q.offset(offset).limit(page_size)
    results = [e.serialize() for e in results]
    result_tuple = (q.count(), results)
    result = create_page(result_tuple, page, page_size, request.url)
    
    Session.close()
    engine.dispose()
    return jsonify(result)

@app.route('/api/v1/<db>/samples')
def samples(db):
    get_db(db)
    db_name = GUDUtils._get_db_name()
    engine, Session = GUDUtils._get_engine_session(db_name)

    table_exists('samples')
    page = check_page(request)
    
    resource = Sample()
    q = resource.select_all_samples(Session)
    offset = (page-1)*page_size
    results = q.offset(offset).limit(page_size)
    results = [e.serialize() for e in results]
    result_tuple = (q.count(), results)
    result = create_page(result_tuple, page, page_size, request.url)
    
    Session.close()
    engine.dispose()
    return jsonify(result)

@app.route('/api/v1/<db>/experiments')
def experiments(db):
    get_db(db)
    db_name = GUDUtils._get_db_name()
    engine, Session = GUDUtils._get_engine_session(db_name)

    table_exists('experiments')
    page = check_page(request)
    
    resource = Experiment()
    q = resource.select_all_experiments(Session)
    offset = (page-1)*page_size
    results = q.offset(offset).limit(page_size)
    results = [e.serialize() for e in results]
    result_tuple = (q.count(), results)
    result = create_page(result_tuple, page, page_size, request.url)
    
    Session.close()
    engine.dispose()
    return jsonify(result)

@app.route('/api/v1/<db>/expression')
def expression(db):
    get_db(db)
    db_name = GUDUtils._get_db_name()
    engine, Session = GUDUtils._get_engine_session(db_name)
    
    table_exists('expression')
    page = check_page(request)

    resource = Expression()
    max_tpm = request.args.get('max_tpm', default=None)    
    min_tpm = request.args.get('min_tpm', default=None)  
    samples = check_split(request.args.get('samples', default=None))  
    q = resource.select_all(Session)

    if max_tpm is not None and min_tpm is not None: 
        try:
            max_tpm = float(max_tpm)
            min_tpm = float(min_tpm)
            q = resource.select_by_expression(Session, min_tpm, max_tpm, q)
        except:
            raise BadRequest('max and min tpm must be floats')
    if samples is not None:
        q = resource.select_by_samples(Session, samples, q)
    
    offset = (page-1)*page_size
    results = q.offset(offset).limit(page_size)
    results = [resource.serialize(e) for e in results]
    result_tuple = (q.count(), results)
    result = create_page(result_tuple, page, page_size, request.url)

    Session.close()
    engine.dispose()
    return jsonify(result)
