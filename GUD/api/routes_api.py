# instructions for adding more API Routes start with '# API_ADDITION(step):'
from GUD.api import app, get_engine_session
from flask import request, jsonify
# API_ADDITION(1): import feature that you would like to add
from GUD.ORM import (Gene, ShortTandemRepeat, CNV, ClinVar, Conservation,
                     DNAAccessibility, Enhancer, HistoneModification, TAD, 
                     TFBinding, TSS, Chrom, Sample, Experiment, Source, Expression) 
from GUD.api.api_helpers import *
from werkzeug.exceptions import BadRequest

# templates 
# API_ADDITION(3): add feature method implementation to query specific feature
# if feature is extention of GF1/2 then follow these templates, else create custom query 
# def GF_feature_template(request, session):                           
#     """retrieves Genomic Feature 2"""
#     resource = Feature()                                                        # instansiate feature 
#     param1 = check_split(request.args.get('param1', default=None))              # get any additional parameters
#     q = genomic_feature_mixin1_queries(session, resource, request)              # build basic query for genomic_feature_mixin1
#     q = genomic_feature_mixin2_queries(session, resource, request, q)           # add to query for genomic_feature_mixin2
#     if param1 is not None:                                                      # for each additional query that is feature specific
#                                                                                 # add to the query                
#         q = resource.select_by_param1(q, param1)
#     return get_result_from_query(q, request, resource, page_size=20, result_tuple_type="genomic_feature") #fetch the query 

# simple resources
def chroms(request, session):                                           
    """retrieves all chromosomes"""
    resource = Chrom()
    q = resource.select_all_chroms(session)
    results = [e.serialize() for e in q]
    json = {}
    if len(results) == 0:
        raise NotFound('No results from this query')
    json = {'results': results}
    return jsonify(json)

def experiments(request, session): 
    """retrieves all experiments"""                                     
    resource = Experiment()
    q = resource.select_all_experiments(session)
    return get_result_from_query(q, request, resource, page_size=20, result_tuple_type="simple")

def samples(request, session):    
    """retrieves all samples"""                                     
    resource = Sample()
    q = resource.select_all_samples(session)
    return get_result_from_query(q, request, resource, page_size=20, result_tuple_type="simple")

def sources(request, session):     
    """retrieves all sources"""                                    
    resource = Source()
    q = resource.select_all_sources(session)
    return get_result_from_query(q, request, resource, page_size=20, result_tuple_type="simple")

#other resources
# TODO: fix this one
# def expression(request, session, page):                                       
#     resource = Expression()
#     max_tpm = request.args.get('max_tpm', default=None)    
#     min_tpm = request.args.get('min_tpm', default=None)  
#     samples = check_split(request.args.get('samples', default=None))  
#     q = resource.select_all(session)
#     if max_tpm is not None and min_tpm is not None: 
#         try:
#             max_tpm = float(max_tpm)
#             min_tpm = float(min_tpm)
#             q = resource.select_by_expression(Session, min_tpm, max_tpm, q)
#         except:
#             raise BadRequest('max and min tpm must be floats')
#     if samples is not None:
#         q = resource.select_by_samples(session, samples, q)
#     page_size = 20
#     offset = (page-1)*page_size
#     results = q.offset(offset).limit(page_size)
#     results = [resource.serialize(e) for e in results]
#     result_tuple = (q.count(), results)
#     result = create_page(result_tuple, page, page_size, request.url)           
#     return jsonify(result)

#GF1 queries 
def clinvar(request, session): 
    """retrieves clinvar variants"""
    resource = ClinVar()
    clinvarIDs = check_split(request.args.get('clinvar_ids', default=None), True)
    q = genomic_feature_mixin1_queries(session, resource, request)
    if clinvarIDs is not None:
        q = resource.select_by_clinvarID(session, q, clinvarIDs)                                               
    return get_result_from_query(q, request, resource, page_size=20, result_tuple_type="genomic_feature")

def conservation(request, session): 
    """retrieves conserved elements"""
    resource = Conservation()
    q = genomic_feature_mixin1_queries(session, resource, request)                                           
    return get_result_from_query(q, request, resource, page_size=20, result_tuple_type="genomic_feature")

def copy_number_variants(request, session):                             
    """retrieves copy number variants"""
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
    return get_result_from_query(q, request, resource, page_size=20, result_tuple_type="genomic_feature")

def genes(request, session):                                           
    """retrieves genes"""
    resource = Gene()
    q = genomic_feature_mixin1_queries(session, resource, request)
    names = check_split(request.args.get('names', default=None))
    if names is not None:
        q = resource.select_by_names(session, q, names)
    return get_result_from_query(q, request, resource, page_size=20, result_tuple_type="genomic_feature")

def short_tandem_repeats(request, session):                            
    """retrieves all STRs"""
    resource = ShortTandemRepeat()
    q = genomic_feature_mixin1_queries(session, resource, request)
    try:
        pathogenic = request.args.get('pathogenic', default=False, type=bool)
        motif = request.args.get('motif', default=None, type=str)
        rotation = request.args.get('rotation', default=False, type=bool)
    except:
        raise BadRequest(
            'rotation must be set to True or False if motif is given')
    if pathogenic:
        q = resource.select_by_pathogenicity(session, q)
    if motif is not None:
        q = resource.select_by_motif(session, motif, q, rotation)
    return get_result_from_query(q, request, resource, page_size=20, result_tuple_type="genomic_feature")

# GF2 queries 
def dna_accessibility(request, session):                                
    resource = DNAAccessibility()
    q = genomic_feature_mixin1_queries(session, resource, request)
    q = genomic_feature_mixin2_queries(session, resource, request, q)
    return get_result_from_query(q, request, resource, page_size=20, result_tuple_type="genomic_feature")

def enhancers(request, session):                                        
    resource = Enhancer()
    q = genomic_feature_mixin1_queries(session, resource, request)
    q = genomic_feature_mixin2_queries(session, resource, request, q)
    return get_result_from_query(q, request, resource, page_size=20, result_tuple_type="genomic_feature")

def histone_modifications(request, session):                           
    resource = HistoneModification()
    histone_types = check_split(request.args.get('histone_types', default=None))
    q = genomic_feature_mixin1_queries(session, resource, request)
    q = genomic_feature_mixin2_queries(session, resource, request, q)
    if histone_types is not None: 
        q = resource.select_by_histone_type(q, histone_types)
    return get_result_from_query(q, request, resource, page_size=20, result_tuple_type="genomic_feature")

def tads(request, session):                                           
    resource = TAD()
    restriction_enzymes = check_split(request.args.get('restriction_enzymes', default=None))
    q = genomic_feature_mixin1_queries(session, resource, request)
    q = genomic_feature_mixin2_queries(session, resource, request, q)
    if restriction_enzymes is not None:
        q = resource.select_by_restriction_enzymes(q, restriction_enzymes)
    return get_result_from_query(q, request, resource, page_size=20, result_tuple_type="genomic_feature")

def tf_binding(request, session):    
    resource = TFBinding()
    tfs = check_split(request.args.get('tfs', default=None))
    q = genomic_feature_mixin1_queries(session, resource, request)
    q = genomic_feature_mixin2_queries(session, resource, request, q)
    if tfs is not None:
        q = resource.select_by_tf(q, tfs)                                
    return get_result_from_query(q, request, resource, page_size=20, result_tuple_type="genomic_feature")

def tss(request, session):
    resource = TSS()
    samples = request.args.get('samples', default=None, type=str)
    if samples is not None:
        return BadRequest('Cannot query TSS table by sample')
    genes = check_split(request.args.get('genes', default=None))
    q = genomic_feature_mixin1_queries(session, resource, request)
    q = genomic_feature_mixin2_queries(session, resource, request, q)
    if genes is not None: 
        q = resource.select_by_genes(q, genes)               
    return get_result_from_query(q, request, resource, page_size=20, result_tuple_type="genomic_feature")

@app.route('/api/v1/<db>/<resource>')
def resource_query(db, resource): 
    """ main control switch function for all valid resources"""
    # API_ADDITION(2): add feature method to switch for querying feature 
    switch = {
        "chroms": chroms,
        "clinvar": clinvar, 
        "copy_number_variants": copy_number_variants,
        "conservation": conservation,
        "dna_accessibility": dna_accessibility,
        "enhancers": enhancers,
        "experiments": experiments,
        # "expression": expression,
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
    response = func(request, Session)   
    Session.close()
    engine.dispose()                              
    return response

# custom control routes

# @app.route('/api/v1/<db>/tss/genic')
# def genic_tss(db):
#     """custom control route getting all genic tss"""
#     engine, Session = get_engine_session(db)
#     table_exists('transcription_start_sites', engine)
#     resource = TSS()
#     q = resource.select_all_genic_tss(Session)   
#     Session.close()
#     engine.dispose()    
#     return get_result_from_query(q, request,  result_tuple_type="genomic_feature", resource=resource)
