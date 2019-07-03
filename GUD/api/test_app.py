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

def test_copy_number_variants():
    #select_by_location
    with app.test_client() as c:
        resp = c.get('/api/v1/copy_number_variants?chrom=chr1&start=852862&end=1008286&location=within')
        data = json.loads(resp.data)
        assert len(data['results']) == 20
        assert data['results'][0]['id'] == 	"copy_number_variants_nssv1609233"
        assert data['size'] == 116
    #select_by_exact_location
    with app.test_client() as c:
        resp = c.get('/api/v1/copy_number_variants?chrom=chr1&start=852863&end=1008286&location=exact')
        data = json.loads(resp.data)
        assert len(data['results']) == 2
        assert data['results'][0]['id'] == "copy_number_variants_nssv1609233"
    #select_by_uids
    with app.test_client() as c:
        resp = c.get('/api/v1/copy_number_variants?uids=nssv1609233')
        data = json.loads(resp.data)
        assert len(data['results']) == 1
        assert data['results'][0]['id'] == "copy_number_variants_nssv1609233"
    #select_by_sources
    with app.test_client() as c:
        resp = c.get('/api/v1/copy_number_variants?sources=ISCA-ClinGen&chrom=chr1&start=852862&end=1008286')
        data = json.loads(resp.data)
        assert data['size'] == 116

def test_gene():
    #select_by_location
    with app.test_client() as c:
        resp = c.get('/api/v1/genes?chrom=chr1&start=11869&end=14362&location=within')
        data = json.loads(resp.data)
        assert len(data['results']) == 3
        assert data['results'][0]['id'] == "genes_1"
        assert data['size'] == 3
    #select_by_exact_location
    with app.test_client() as c:
        resp = c.get('/api/v1/genes?chrom=chr1&start=11869&end=14362&location=exact')
        data = json.loads(resp.data)
        assert len(data['results']) == 1
        assert data['results'][0]['id'] == "genes_1"
    #select_by_uids
    with app.test_client() as c:
        resp = c.get('/api/v1/genes?uids=1')
        data = json.loads(resp.data)
        assert len(data['results']) == 1
        assert data['results'][0]['id'] == "genes_1"
    #select_by_sources
    with app.test_client() as c:
        resp = c.get('/api/v1/genes?chrom=chr1&start=11869&end=14362&sources=refGene')
        data = json.loads(resp.data)
        assert data['size'] == 3
    #select_by_names
    with app.test_client() as c:
        resp = c.get('/api/v1/genes?names=DDX11L1')
        data = json.loads(resp.data)
        assert data['size'] == 1
    #get_all_gene_symbols
    with app.test_client() as c:
        resp = c.get('/api/v1/genesymbols')
        data = json.loads(resp.data)
        assert data['size'] == 28194

def test_short_tandem_repeat():
#select_by_location
    with app.test_client() as c:
        resp = c.get('/api/v1/short_tandem_repeats?chrom=chr1&start=11869&end=14362&location=within')
        data = json.loads(resp.data)
        assert len(data['results']) == 1
        assert data['results'][0]['id'] == "short_tandem_repeats_305056"
        assert data['size'] == 1
    #select_by_exact_location
    with app.test_client() as c:
        resp = c.get('/api/v1/short_tandem_repeats?chrom=chr1&start=14070&end=14081&location=within')
        data = json.loads(resp.data)
        assert len(data['results']) == 1
        assert data['results'][0]['id'] == "short_tandem_repeats_305056"
    #select_by_uids
    with app.test_client() as c:
        resp = c.get('/api/v1/short_tandem_repeats?uids=305056')
        data = json.loads(resp.data)
        assert len(data['results']) == 1
        assert data['results'][0]['id'] == "short_tandem_repeats_305056"
    #select_by_sources
    with app.test_client() as c:
        resp = c.get('/api/v1/short_tandem_repeats?chrom=chr1&start=11869&end=14362&sources=gangSTR')
        data = json.loads(resp.data)
        assert data['size'] == 1
    #select_by_pathogencity 
    with app.test_client() as c:
        resp = c.get('/api/v1/short_tandem_repeats?pathogenicity=True')
        data = json.loads(resp.data)
        assert data['size'] == 26
    #select_by_motif
    with app.test_client() as c:
        resp = c.get('/api/v1/short_tandem_repeats?motif=ATGGG&rotation=True')
        data = json.loads(resp.data)
        assert data['size'] == 92
    

# def test_conservation(client):
#     #select_by_location
#     #select_by_exact_location
#     #select_by_uids
#     #select_by_sources
#     return False

# def test_repeat_mask():
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
