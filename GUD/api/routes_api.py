from GUD import GUDUtils
from GUD.api import app
from flask import request, jsonify
from GUD.ORM import (Gene, ShortTandemRepeat, CNV, ClinVar, Conservation,
                     DNAAccessibility, Enhancer, HistoneModification, TAD, 
                     TFBinding, TSS, Chrom, Sample, Experiment, Source, Expression)
from werkzeug.exceptions import HTTPException, NotFound, BadRequest
import sys
from GUD.api.api_helpers import *
page_size = 20
# print(names, file=sys.stdout)

@app.route('/api/v1/clinvar')
def clinvar():
    page = check_page(request)
    session = GUDUtils.get_session()
    resource = ClinVar()
    clinvarIDs = check_split(request.args.get(
        'clinvar_ids', default=None), True)
    q = genomic_feature_mixin1_queries(session, resource, request)
    if clinvarIDs is not None:
        q = resource.select_by_clinvarID(session, q, clinvarIDs)
    
    session.close()
    result_tuple = get_genomic_feature_results(resource, q, page)
    result = create_page(result_tuple, page, request.url)
    return jsonify(result)


@app.route('/api/v1/copy_number_variants')
def copy_number_variants():
    page = check_page(request)
    
    session = GUDUtils.get_session()
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
    session.close()
    result_tuple = get_genomic_feature_results(resource, q, page)
    result = create_page(result_tuple, page, request.url)
    return jsonify(result)


@app.route('/api/v1/genes')
def genes():
    page = check_page(request)
    
    session = GUDUtils.get_session()
    resource = Gene()
    q = genomic_feature_mixin1_queries(session, resource, request)
    names = check_split(request.args.get('names', default=None))
    if names is not None:
        q = resource.select_by_names(session, q, names)
    
    session.close()
    result_tuple = get_genomic_feature_results(resource, q, page)
    result = create_page(result_tuple, page, request.url)
    return jsonify(result)


@app.route('/api/v1/genes/symbols')
def gene_symbols():
    url = request.url
    
    session = GUDUtils.get_session()
    page = int(request.args.get('page', default=1))
    q = Gene().get_all_gene_symbols(session)
    
    session.close()
    offset = (page-1)*page_size
    result = q.offset(offset).limit(page_size)
    result =  [g[0] for g in result]
    result_tuple = (q.count(), result)
    result = create_page(result_tuple, page, request.url)                                 
    return jsonify(result)


@app.route('/api/v1/short_tandem_repeats')
def strs():
    page = check_page(request)
    
    session = GUDUtils.get_session()
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
    result_tuple = get_genomic_feature_results(resource, q, page)
    result = create_page(result_tuple, page, request.url)
    
    session.close()
    return jsonify(result)


@app.route('/api/v1/short_tandem_repeats/pathogenic')
def pathogenic_strs():
    page = check_page(request)
    
    session = GUDUtils.get_session()
    resource = ShortTandemRepeat()
    q = resource.select_by_pathogenicity(session)
    result_tuple = get_genomic_feature_results(resource, q, page)
    result = create_page(result_tuple, page, request.url)
    
    session.close()
    return jsonify(result)

@app.route('/api/v1/enhancers')
@app.route('/api/v1/dna_accessibility')
def mixin2():
    page = check_page(request)
    
    session = GUDUtils.get_session()
    if (request.path == '/api/v1/dna_accessibility'):
        resource = DNAAccessibility()
    elif (request.path == '/api/v1/enhancers'):
        resource = Enhancer()
    q = genomic_feature_mixin1_queries(session, resource, request)
    q = genomic_feature_mixin2_queries(session, resource, request, q)
    
    session.close()
    result_tuple = get_genomic_feature_results(resource, q, page)
    result = create_page(result_tuple, page, request.url)
    return jsonify(result)

@app.route('/api/v1/histone_modifications')
def histone_modifications():
    page = check_page(request)
    
    session = GUDUtils.get_session()
    resource = HistoneModification()
    histone_types = check_split(request.args.get('histone_types', default=None))
    q = genomic_feature_mixin1_queries(session, resource, request)
    q = genomic_feature_mixin2_queries(session, resource, request, q)
    if histone_types is not None: 
        q = resource.select_by_histone_type(q, histone_types)
    
    session.close()
    result_tuple = get_genomic_feature_results(resource, q, page)
    result = create_page(result_tuple, page, request.url)
    return jsonify(result)

@app.route('/api/v1/tads')
def tads():
    page = check_page(request)
    
    session = GUDUtils.get_session()
    resource = TAD()
    restriction_enzymes = check_split(request.args.get('restriction_enzymes', default=None))
    q = genomic_feature_mixin1_queries(session, resource, request)
    q = genomic_feature_mixin2_queries(session, resource, request, q)
    if restriction_enzymes is not None:
        q = resource.select_by_restriction_enzymes(q, restriction_enzymes)
    
    session.close()
    result_tuple = get_genomic_feature_results(resource, q, page)
    result = create_page(result_tuple, page, request.url)
    return jsonify(result)

@app.route('/api/v1/tf_binding')
def tf_binding():
    page = check_page(request)
    
    session = GUDUtils.get_session()
    resource = TFBinding()
    tfs = check_split(request.args.get('tfs', default=None))
    q = genomic_feature_mixin1_queries(session, resource, request)
    q = genomic_feature_mixin2_queries(session, resource, request, q)
    if tfs is not None:
        q = resource.select_by_tf(q, tfs)
    
    session.close()
    result_tuple = get_genomic_feature_results(resource, q, page)
    result = create_page(result_tuple, page, request.url)
    return jsonify(result)

@app.route('/api/v1/tss')
def tss():
    samples = request.args.get('samples', default=None, type=str)
    if samples is not None:
        return BadRequest('Cannot query TSS table by sample')
    page = check_page(request)
    
    session = GUDUtils.get_session()
    resource = TSS()
    genes = check_split(request.args.get('genes', default=None))
    q = genomic_feature_mixin1_queries(session, resource, request)
    q = genomic_feature_mixin2_queries(session, resource, request, q)
    if genes is not None: 
        q = resource.select_by_genes(q, genes)
    
    session.close()
    result_tuple = get_genomic_feature_results(resource, q, page)
    result = create_page(result_tuple, page, request.url)
    return jsonify(result)

@app.route('/api/v1/tss/genic')
def genic_tss():
    page = check_page(request)
    
    session = GUDUtils.get_session()
    resource = TSS()
    q = resource.select_all_genic_tss(session)
    result_tuple = get_genomic_feature_results(resource, q, page)
    result = create_page(result_tuple, page, request.url)
    return jsonify(result)

@app.route('/api/v1/chroms')
def chroms():
    page = check_page(request)
    
    session = GUDUtils.get_session()
    resource = Chrom()
    q = resource.select_all_chroms(session)
    offset = (page-1)*page_size
    results = q.offset(offset).limit(page_size)
    results = [e.serialize() for e in results]
    result_tuple = (q.count(), results)
    result = create_page(result_tuple, page, request.url)
    
    session.close()
    return jsonify(result)

@app.route('/api/v1/sources')
def sources():
    page = check_page(request)
    
    session = GUDUtils.get_session()
    resource = Source()
    q = resource.select_all_sources(session)
    offset = (page-1)*page_size
    results = q.offset(offset).limit(page_size)
    results = [e.serialize() for e in results]
    result_tuple = (q.count(), results)
    result = create_page(result_tuple, page, request.url)
    
    session.close()
    return jsonify(result)

@app.route('/api/v1/samples')
def samples():
    page = check_page(request)
    
    session = GUDUtils.get_session()
    resource = Sample()
    q = resource.select_all_samples(session)
    offset = (page-1)*page_size
    results = q.offset(offset).limit(page_size)
    results = [e.serialize() for e in results]
    result_tuple = (q.count(), results)
    result = create_page(result_tuple, page, request.url)
    
    session.close()
    return jsonify(result)

@app.route('/api/v1/experiments')
def experiments():
    page = check_page(request)
    
    session = GUDUtils.get_session()
    resource = Experiment()
    q = resource.select_all_experiments(session)
    offset = (page-1)*page_size
    results = q.offset(offset).limit(page_size)
    results = [e.serialize() for e in results]
    result_tuple = (q.count(), results)
    result = create_page(result_tuple, page, request.url)
    
    session.close()
    return jsonify(result)

@app.route('/api/v1/expression')
def expression():
    page = check_page(request)
    
    session = GUDUtils.get_session()
    resource = Expression()
    q = None
    max_tpm = request.args.get('max_tpm', default=None)    
    min_tpm = request.args.get('min_tpm', default=None)  
    samples = check_split(request.args.get('samples', default=None))  
    if max_tpm is not None and min_tpm is not None: 
        try:
            max_tpm = float(max_tpm)
            min_tpm = float(min_tpm)
            q = resource.select_by_expression(session, min_tpm, max_tpm, q)
        except:
            raise BadRequest('max and min tpm must be floats')
    if samples is not None:
        q = resource.select_by_samples(session, samples, q)

    offset = (page-1)*page_size
    results = q.offset(offset).limit(page_size)
    results = [resource.serialize(e) for e in results]
    result_tuple = (q.count(), results)
    result = create_page(result_tuple, page, request.url)
    
    session.close()
    return jsonify(result)
