from GUD import GUDUtils
from GUD.api import app, get_engine_session
from flask import request, jsonify
from GUD.ORM import (Gene, ShortTandemRepeat, CNV, ClinVar, Conservation,
                     DNAAccessibility, Enhancer, HistoneModification, TAD, 
                     TFBinding, TSS, Chrom, Sample, Experiment, Source, Expression)
from GUD.api.api_helpers import *
from werkzeug.exceptions import HTTPException, NotFound, BadRequest
import sys
# page_size = 20

## move to api helpers 
def get_result_tuple_simple(query, page, page_size=20):
    offset = (page-1)*page_size
    results = query.offset(offset).limit(page_size)
    results = [e.serialize() for e in results]
    return(query.count(), results)

def get_result_from_query(query, page, page_size=20, result_tuple_type="simple", resource=None):
    if result_tuple_type == "simple": 
        result_tuple = get_result_tuple_simple(query, page)
    elif result_tuple_type == "genomic_feature": 
        result_tuple = get_genomic_feature_results(resource, query, page_size,  page)
    result = create_page(result_tuple, page, page_size, request.url)           
    return jsonify(result)

# simple resources
def chroms(request, session, page):                                           
    resource = Chrom()
    q = resource.select_all_chroms(session)
    return get_result_from_query(q, page)

def experiments(request, session, page):                                      
    resource = Experiment()
    q = resource.select_all_experiments(session)
    return get_result_from_query(q, page)

def samples(request, session, page):                                         
    resource = Sample()
    q = resource.select_all_samples(session)
    return get_result_from_query(q, page)

def sources(request, session, page):                                         
    resource = Source()
    q = resource.select_all_sources(session)
    return get_result_from_query(q, page)

#other resources
def expression(request, session, page):                                       
    return False

#GF1 queries 
def clinvar(request, session, page): 
    resource = ClinVar()
    clinvarIDs = check_split(request.args.get('clinvar_ids', default=None), True)
    q = genomic_feature_mixin1_queries(session, resource, request)
    if clinvarIDs is not None:
        q = resource.select_by_clinvarID(session, q, clinvarIDs)                                               
    return get_result_from_query(q, page, result_tuple_type="genomic_feature", resource=resource)

def copy_number_variants(request, session, page):                             
    clinical_assertion  = request.args.get('clinical_assertion', default=None)
    clinvar_accession   = request.args.get('clinvar_accession', default=None)
    dbVar_accession     = request.args.get('dbVar_accession', default=None)
    resource = CNV()
    q = genomic_feature_mixin1_queries(session, resource, request)
    if clinical_assertion is not None:
        q = resource.select_by_clinical_assertion(session, q, clinical_assertion)
    if clinvar_accession is not None:
        q = resource.select_by_clinvar_accession(session, q, clinvar_accession)
    if dbVar_accession is not None:
        q = resource.select_by_dbvar_accession(session, q, dbVar_accession)
    return get_result_from_query(q, page, result_tuple_type="genomic_feature", resource=resource)

def genes(request, session, page):                                           
    resource = Gene()
    q = genomic_feature_mixin1_queries(session, resource, request)
    names = check_split(request.args.get('names', default=None))
    if names is not None:
        q = resource.select_by_names(session, q, names)
    if (q is None):
        raise BadRequest('Query not specified correctly')
    return get_result_from_query(q, page, result_tuple_type="genomic_feature", resource=resource)

def short_tandem_repeats(request, session, page):                            
    resource = ShortTandemRepeat()
    q = genomic_feature_mixin1_queries(session, resource, request)
    try:
        motif = request.args.get('motif', default=None, type=str)
        rotation = request.args.get('rotation', default=False, type=bool)
    except:
        raise BadRequest(
            'rotation must be set to True or False if motif is given')
    if motif is not None:
        q = resource.select_by_motif(session, motif, q, rotation)
    return get_result_from_query(q, page, result_tuple_type="genomic_feature", resource=resource)

# GF2 queries 
def dna_accessibility(request, session, page):                                
    return False

def enhancers(request, session, page):                                        
    return False

def histone_modifications(request, session, page):                           
    return False

def tads(request, session, page):                                           
    return False

def tf_binding(request, session, page):                                     
    return False

def tss(request, session, page):                                           
    return False

@app.route('/api/v1/<db>/<resource>')
def resource_query(db, resource): 
    switch = {
        "chroms": chroms,
        "clinvar": clinvar, 
        "copy_number_variants": copy_number_variants,
        "dna_accessibility": dna_accessibility,
        "enhancers": enhancers,
        "experiments": experiments,
        "expression": expression,
        "genes": genes,
        "histone_modifications": histone_modifications,
        "samples": samples,
        "short_tandem_repeats": short_tandem_repeats,
        "sources": sources,
        "tads": tads,
        "tf_binding": tf_binding,
        "tss": tss
    }
    engine, Session = get_engine_session(db)
    func = switch.get(resource, "none")
    if func == "none":                  # check if this is invalid route 
        raise BadRequest('Invalid resource')
    
    table_exists(resource, engine)      # check that table exists 
    page = check_page(request)          # check that page is valid 
    response = func(request, Session, page)                                    
    return response



# @app.route('/api/v1/<db>/genes/symbols')
# def gene_symbols(db):
#     Session = get_session(db)

#     table_exists('genes', get_engine(db))
#     url = request.url
#     page = int(request.args.get('page', default=1))
#     q = Gene().get_all_gene_symbols(Session)
    
#     page_size = 1000 ## set custom page size 
#     offset = (page-1)*page_size
#     result = q.offset(offset).limit(page_size)
#     result =  [g[0] for g in result]
    
#     result_tuple = (q.count(), result)
#     result = create_page(result_tuple, page, page_size, request.url)           
#     page_size = 20  ## reset page size
#     response = jsonify(result)                                                 
#     return response

# @app.route('/api/v1/<db>/short_tandem_repeats/pathogenic')
# def pathogenic_strs(db):
#     Session = get_session(db)

#     table_exists('short_tandem_repeats', get_engine(db))
#     page = check_page(request)
    
#     resource = ShortTandemRepeat()
#     q = resource.select_by_pathogenicity(Session)
    
#     result_tuple = get_genomic_feature_results(resource, q, page_size,  page)
#     result = create_page(result_tuple, page, page_size, request.url)           #same
#     return jsonify(result)

# @app.route('/api/v1/<db>/enhancers')
# @app.route('/api/v1/<db>/dna_accessibility')
# def mixin2(db):
#     Session = get_session(db)
  
#     page = check_page(request)
    
#     if ('dna_accessibility' in request.path):
#         table_exists('dna_accessibility', get_engine(db))
#         resource = DNAAccessibility()
#     elif ('enhancers' in request.path):
#         table_exists('enhancers', get_engine(db))
#         resource = Enhancer()
#     q = genomic_feature_mixin1_queries(Session, resource, request)
#     q = genomic_feature_mixin2_queries(Session, resource, request, q)

#     result_tuple = get_genomic_feature_results(resource, q, page_size,  page)
#     result = create_page(result_tuple, page, page_size, request.url)           #same
#     return jsonify(result)

# @app.route('/api/v1/<db>/histone_modifications')
# def histone_modifications(db):
#     Session = get_session(db)
 
#     table_exists('histone_modifications', get_engine(db))
#     page = check_page(request)

#     resource = HistoneModification()
#     histone_types = check_split(request.args.get('histone_types', default=None))
#     q = genomic_feature_mixin1_queries(Session, resource, request)
#     q = genomic_feature_mixin2_queries(Session, resource, request, q)

#     if histone_types is not None: 
#         q = resource.select_by_histone_type(q, histone_types)

#     result_tuple = get_genomic_feature_results(resource, q, page_size,  page)
#     result = create_page(result_tuple, page,page_size, request.url)                                 
#     return jsonify(result)

# @app.route('/api/v1/<db>/tads')
# def tads(db):
#     Session = get_session(db)
#     table_exists('tads', get_engine(db))
#     page = check_page(request)
    
#     resource = TAD()
#     restriction_enzymes = check_split(request.args.get('restriction_enzymes', default=None))
#     q = genomic_feature_mixin1_queries(Session, resource, request)
#     q = genomic_feature_mixin2_queries(Session, resource, request, q)
#     if restriction_enzymes is not None:
#         q = resource.select_by_restriction_enzymes(q, restriction_enzymes)

#     result_tuple = get_genomic_feature_results(resource, q, page_size,  page)
#     result = create_page(result_tuple, page, page_size, request.url)           
#     return jsonify(result)

# @app.route('/api/v1/<db>/tf_binding')
# def tf_binding(db):
#     Session = get_session(db)

#     table_exists('tf_binding', get_engine(db))
#     page = check_page(request)
    
#     resource = TFBinding()
#     tfs = check_split(request.args.get('tfs', default=None))
#     q = genomic_feature_mixin1_queries(Session, resource, request)
#     q = genomic_feature_mixin2_queries(Session, resource, request, q)
#     if tfs is not None:
#         q = resource.select_by_tf(q, tfs)
    
#     result_tuple = get_genomic_feature_results(resource, q, page_size,  page)
#     result = create_page(result_tuple, page, page_size, request.url)           
#     return jsonify(result)

# @app.route('/api/v1/<db>/tss')
# def tss(db):
#     Session = get_session(db)

#     table_exists('transcription_start_sites', get_engine(db))
#     samples = request.args.get('samples', default=None, type=str)
#     if samples is not None:
#         return BadRequest('Cannot query TSS table by sample')
#     page = check_page(request)
    
#     resource = TSS()
#     genes = check_split(request.args.get('genes', default=None))
#     q = genomic_feature_mixin1_queries(Session, resource, request)
#     q = genomic_feature_mixin2_queries(Session, resource, request, q)
#     if genes is not None: 
#         q = resource.select_by_genes(q, genes)

#     result_tuple = get_genomic_feature_results(resource, q, page_size,  page)
#     result = create_page(result_tuple, page, page_size, request.url)           
#     return jsonify(result)

# @app.route('/api/v1/<db>/tss/genic')
# def genic_tss(db):
#     Session = get_session(db)
#     table_exists('transcription_start_sites', get_engine(db))
#     page = check_page(request)
#     resource = TSS()
#     q = resource.select_all_genic_tss(Session)
    
#     result_tuple = get_genomic_feature_results(resource, q, page_size,  page)
#     result = create_page(result_tuple, page, page_size, request.url)           
#     return jsonify(result)

# @app.route('/api/v1/<db>/expression')
# def expression(db):
#     Session = get_session(db)
    
#     table_exists('expression', get_engine(db))
#     page = check_page(request)

#     resource = Expression()
#     max_tpm = request.args.get('max_tpm', default=None)    
#     min_tpm = request.args.get('min_tpm', default=None)  
#     samples = check_split(request.args.get('samples', default=None))  
#     q = resource.select_all(Session)

#     if max_tpm is not None and min_tpm is not None: 
#         try:
#             max_tpm = float(max_tpm)
#             min_tpm = float(min_tpm)
#             q = resource.select_by_expression(Session, min_tpm, max_tpm, q)
#         except:
#             raise BadRequest('max and min tpm must be floats')
#     if samples is not None:
#         q = resource.select_by_samples(Session, samples, q)
    
#     offset = (page-1)*page_size
#     results = q.offset(offset).limit(page_size)
#     results = [resource.serialize(e) for e in results]
#     result_tuple = (q.count(), results)
#     result = create_page(result_tuple, page, page_size, request.url)           
#     return jsonify(result)
