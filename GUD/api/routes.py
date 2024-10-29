from GUD.api import app, get_engine_session
from flask import request, jsonify
from GUD.ORM import (Gene, ShortTandemRepeat, CNV, ClinVar, Chrom) 
from GUD.api.api_helpers import *
from werkzeug.exceptions import BadRequest
from GUD.api import app
from flask import request, jsonify
from werkzeug.exceptions import NotFound, BadRequest
import json,os, sys

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

def sources(request, session):     
    """retrieves all sources"""                                    
    resource = Source()
    q = resource.select_all_sources(session)
    return get_result_from_query(q, request, resource, result_tuple_type="simple")

#GF1 queries 
def clinvar(request, session): 
    """retrieves clinvar variants"""
    resource = ClinVar()
    clinvarIDs = check_split(request.args.get('clinvar_ids', default=None), True)
    q = genomic_feature_mixin1_queries(session, resource, request)
    if clinvarIDs is not None:
        q = [resource.select_by_clinvarID(session, i, clinvarIDs) for i in q]                                            
    return get_result_from_query(q, request, resource, result_tuple_type="genomic_feature")

def copy_number_variants(request, session):                             
    """retrieves copy number variants"""
    resource = CNV()
    q = genomic_feature_mixin1_queries(session, resource, request)
    return get_result_from_query(q, request, resource, result_tuple_type="genomic_feature")

def genes(request, session):                                           
    """retrieves genes"""
    resource = Gene()
    q= genomic_feature_mixin1_queries(session, resource, request)
    names = check_split(request.args.get('names', default=None))
    if names is not None:
        q = [resource.select_by_names(session, i, names) for i in q]  
    return get_result_from_query(q, request, resource,
                                result_tuple_type="genomic_feature")

def short_tandem_repeats(request, session):                            
    """retrieves all STRs"""
    resource = ShortTandemRepeat()
    q= genomic_feature_mixin1_queries(session, resource, request)
    return get_result_from_query(q, request, resource, 
                                result_tuple_type="genomic_feature")

import time
@app.route('/api/v1/<db>/<resource>')
def resource_query(db, resource): 
    switch = {
    "chroms": chroms,
    "clinvar": clinvar, 
    "copy_number_variants": copy_number_variants,
    "genes": genes,
    "short_tandem_repeats": short_tandem_repeats,
    "sources": sources,
    }
    engine, session = get_engine_session(db)
    func = switch.get(resource, "none")
    if func == "none":                  # check if this is invalid route 
        raise BadRequest('Invalid resource')
    table_exists(resource, engine)      # check that table exists 
    response = func(request, session)
    session.close()
    engine.dispose()
    return response