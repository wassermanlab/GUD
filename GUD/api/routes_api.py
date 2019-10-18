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
    resource = Expression()
    max_tpm = request.args.get('max_tpm', default=None)    
    min_tpm = request.args.get('min_tpm', default=None)  
    samples = check_split(request.args.get('samples', default=None))  
    q = resource.select_all(session)
    if max_tpm is not None and min_tpm is not None: 
        try:
            max_tpm = float(max_tpm)
            min_tpm = float(min_tpm)
            q = resource.select_by_expression(Session, min_tpm, max_tpm, q)
        except:
            raise BadRequest('max and min tpm must be floats')
    if samples is not None:
        q = resource.select_by_samples(session, samples, q)
    page_size = 20
    offset = (page-1)*page_size
    results = q.offset(offset).limit(page_size)
    results = [resource.serialize(e) for e in results]
    result_tuple = (q.count(), results)
    result = create_page(result_tuple, page, page_size, request.url)           
    return jsonify(result)

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
    resource = DNAAccessibility()
    q = genomic_feature_mixin1_queries(session, resource, request)
    q = genomic_feature_mixin2_queries(session, resource, request, q)
    return get_result_from_query(q, page, result_tuple_type="genomic_feature", resource=resource)

def enhancers(request, session, page):                                        
    resource = Enhancer()
    q = genomic_feature_mixin1_queries(session, resource, request)
    q = genomic_feature_mixin2_queries(session, resource, request, q)
    return get_result_from_query(q, page, result_tuple_type="genomic_feature", resource=resource)

def histone_modifications(request, session, page):                           
    resource = HistoneModification()
    histone_types = check_split(request.args.get('histone_types', default=None))
    q = genomic_feature_mixin1_queries(session, resource, request)
    q = genomic_feature_mixin2_queries(session, resource, request, q)
    if histone_types is not None: 
        q = resource.select_by_histone_type(q, histone_types)
    return get_result_from_query(q, page, result_tuple_type="genomic_feature", resource=resource)

def tads(request, session, page):                                           
    resource = TAD()
    restriction_enzymes = check_split(request.args.get('restriction_enzymes', default=None))
    q = genomic_feature_mixin1_queries(session, resource, request)
    q = genomic_feature_mixin2_queries(session, resource, request, q)
    if restriction_enzymes is not None:
        q = resource.select_by_restriction_enzymes(q, restriction_enzymes)
    return get_result_from_query(q, page, result_tuple_type="genomic_feature", resource=resource)

def tf_binding(request, session, page):    
    resource = TFBinding()
    tfs = check_split(request.args.get('tfs', default=None))
    q = genomic_feature_mixin1_queries(session, resource, request)
    q = genomic_feature_mixin2_queries(session, resource, request, q)
    if tfs is not None:
        q = resource.select_by_tf(q, tfs)                                
    return get_result_from_query(q, page, result_tuple_type="genomic_feature", resource=resource)

def tss(request, session, page):
    resource = TSS()
    samples = request.args.get('samples', default=None, type=str)
    if samples is not None:
        return BadRequest('Cannot query TSS table by sample')
    genes = check_split(request.args.get('genes', default=None))
    q = genomic_feature_mixin1_queries(session, resource, request)
    q = genomic_feature_mixin2_queries(session, resource, request, q)
    if genes is not None: 
        q = resource.select_by_genes(q, genes)               
    return get_result_from_query(q, page, result_tuple_type="genomic_feature", resource=resource)

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

@app.route('/api/v1/<db>/genes/symbols')
def gene_symbols(db):
    engine, Session = get_engine_session(db)
    table_exists('genes', engine)
    url = request.url
    page = check_page(request)
    q = Gene().get_all_gene_symbols(Session)
    page_size = 1000
    offset = (page-1)*page_size
    result = q.offset(offset).limit(page_size)
    result =  [g[0] for g in result]
    result_tuple = (q.count(), result)
    result = create_page(result_tuple, page, page_size, request.url)                                            
    return jsonify(result) 

@app.route('/api/v1/<db>/short_tandem_repeats/pathogenic')
def pathogenic_strs(db):
    engine, Session = get_engine_session(db)
    table_exists('short_tandem_repeats', engine)
    page = check_page(request)
    resource = ShortTandemRepeat()
    q = resource.select_by_pathogenicity(Session)
    return get_result_from_query(q, page, result_tuple_type="genomic_feature", resource=resource)

@app.route('/api/v1/<db>/tss/genic')
def genic_tss(db):
    engine, Session = get_engine_session(db)
    table_exists('transcription_start_sites', engine)
    page = check_page(request)
    resource = TSS()
    q = resource.select_all_genic_tss(Session)        
    return get_result_from_query(q, page, result_tuple_type="genomic_feature", resource=resource)
