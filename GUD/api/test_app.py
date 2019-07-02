import os, sys
import tempfile

import pytest

from GUD.api import app
import flask
from flask import request, jsonify, json

# ============== gfmixin1 ============== #

with app.test_request_context('/api/v1/clinvar?clinvar_id=475283'):
    assert flask.request.args['clinvar_id'] == '475283'

def test_clinvar():
    #select_by_location
    with app.test_client() as c:
        resp = c.get('/api/v1/clinvar?location=within&chrom=chr1&start=979722&end=979780')
        data = json.loads(resp.data)
        assert len(data['results']) == 2
        assert data['results'][0]['id'] == "clinvar_105"
        assert data['results'][1]['id'] == "clinvar_106"
    #select_by_exact_location
    with app.test_client() as c:
        resp = c.get('/api/v1/clinvar?chrom=chr1&start=949422&end=949422&location=exact')
        data = json.loads(resp.data)
        assert len(data['results']) == 1
        assert data['results'][0]['id'] == "clinvar_1"
    #select_by_uids
    with app.test_client() as c:
        resp = c.get('/api/v1/clinvar?uids=104')
        data = json.loads(resp.data)
        assert len(data['results']) == 1
        assert data['results'][0]['id'] == "clinvar_104"
    #select_by_sources
    with app.test_client() as c:
        resp = c.get('/api/v1/clinvar?sources=ClinVar_2018-10-28&chrom=chr1&start=979722&end=979780')
        data = json.loads(resp.data)
        assert len(data['results']) == 2
        assert data['results'][0]['id'] == "clinvar_105"
        assert data['results'][1]['id'] == "clinvar_106"
    #select_by_clinvarID
    with app.test_client() as c:
        resp = c.get('/api/v1/clinvar?clinvar_id=475283')
        data = json.loads(resp.data)
        assert len(data['results']) == 1
        assert data['results'][0]['id'] == "clinvar_1"

# def test_conservation(client):
#     #select_by_location
#     #select_by_exact_location
#     #select_by_uids
#     #select_by_sources
#     return False

# def test_copy_number_variants(client):
#     #select_by_location
#     #select_by_exact_location
#     #select_by_uids
#     #select_by_sources
#     return False

# def test_gene(client):
#     #select_by_location
#     #select_by_exact_location
#     #select_by_uids
#     #select_by_sources
#     #get_all_gene_symbols
#     #select_by_names
#     return False

# def test_repeat_mask(client):
#     #select_by_location
#     #select_by_exact_location
#     #select_by_uids
#     #select_by_sources
#     return False

# def test_short_tandem_repeat(client):
#     #select_by_location
#     #select_by_exact_location
#     #select_by_uids
#     #select_by_sources
#     return False

# # ============== gfmixin2 ============== #
# def test_dna_accessibility(client):
#     return False

# def test_enhancer(client):
#     return False

# def test_histone_modification(client):
#     return False

# def test_tad(client):
#     return False

# def test_tf_binding(client):
#     return False

# def test_tss(client):
#     return False

# # ============== other ============== #
# def test_sample(client):
#     return False

# def test_experiment(client):
#     return False

# def test_source(client):
#     return False

# def test_chrom(client):
#     return False
