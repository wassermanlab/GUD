import os, sys
import tempfile

import pytest

from GUD.api import app
import flask
from flask import request, jsonify, json

# ============== gfmixin1 ============== #
@pytest.mark.skip(reason="tamar_test")
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
        resp = c.get('/api/v1/clinvar?uids=105&chrom=chr1&start=979722&end=979780&location=within')##
        data = json.loads(resp.data)
        assert len(data['results']) == 1
        assert data['results'][0]['id'] == "clinvar_105"
    #select_by_sources
    with app.test_client() as c:
        resp = c.get('/api/v1/clinvar?sources=ClinVar_2018-10-28&chrom=chr1&start=979722&end=979780&location=within')
        data = json.loads(resp.data)
        assert len(data['results']) == 2
        assert data['results'][0]['id'] == "clinvar_105"
        assert data['results'][1]['id'] == "clinvar_106"
    #select_by_clinvarID
    with app.test_client() as c:
        resp = c.get('/api/v1/clinvar?clinvar_ids=541159&chrom=chr1&start=979722&end=979780&location=within')##
        data = json.loads(resp.data)
        assert len(data['results']) == 1
        assert data['results'][0]['id'] == "clinvar_105"

@pytest.mark.skip(reason="tamar_test")
def test_copy_number_variants():
    #select_by_location
    with app.test_client() as c:
        resp = c.get('/api/v1/copy_number_variants?chrom=chr1&start=852862&end=1008286&location=overlapping')
        data = json.loads(resp.data)
        assert len(data['results']) == 20
        assert data['results'][0]['id'] == 	"copy_number_variants_nssv13653606"
        assert data['size'] == 116
    #select_by_exact_location
    with app.test_client() as c:
        resp = c.get('/api/v1/copy_number_variants?chrom=chr1&start=852863&end=1008286&location=exact')
        data = json.loads(resp.data)
        assert len(data['results']) == 2
        assert data['results'][0]['id'] == "copy_number_variants_nssv1609233"
    #select_by_uids
    with app.test_client() as c:
        resp = c.get('/api/v1/copy_number_variants?uids=nssv1609233&chrom=chr1&start=852862&end=1008286&location=overlapping')
        data = json.loads(resp.data)
        assert len(data['results']) == 1
        assert data['results'][0]['id'] == "copy_number_variants_nssv1609233"
    #select_by_sources
    with app.test_client() as c:
        resp = c.get('/api/v1/copy_number_variants?sources=ISCA-ClinGen&chrom=chr1&start=852862&end=1008286&location=overlapping')
        data = json.loads(resp.data)
        assert data['size'] == 116

# @pytest.mark.skip(reason="tamar_test")
def test_gene():
    #select_by_location
    with app.test_client() as c:
        resp = c.get('/api/v1/genes?chrom=chr1&start=11869&end=14362&location=overlapping')
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
        resp = c.get('/api/v1/genes?chrom=chr1&start=11869&end=14362&location=exact&uids=1')
        data = json.loads(resp.data)
        assert len(data['results']) == 1
        assert data['results'][0]['id'] == "genes_1"
    #select_by_sources
    with app.test_client() as c:
        resp = c.get('/api/v1/genes?chrom=chr1&start=11869&end=14362&sources=refGene&location=overlapping')
        data = json.loads(resp.data)
        assert data['size'] == 3
    #select_by_names
    with app.test_client() as c:
        resp = c.get('/api/v1/genes?chrom=chr1&start=11869&end=14362&sources=refGene&location=overlapping&names=DDX11L1')
        data = json.loads(resp.data)
        assert data['size'] == 1
    #get_all_gene_symbols
    with app.test_client() as c:
        resp = c.get('/api/v1/genes/symbols')
        data = json.loads(resp.data)
        # assert data['size'] == 28194

@pytest.mark.skip(reason="tamar_test")
def test_short_tandem_repeat():
    #select_by_location
    with app.test_client() as c:
        resp = c.get('/api/v1/short_tandem_repeats?chrom=chr1&start=11869&end=14362&location=overlapping')
        data = json.loads(resp.data)
        assert len(data['results']) == 1
        assert data['results'][0]['id'] == "short_tandem_repeats_305056"
        assert data['size'] == 1
    #select_by_exact_location
    with app.test_client() as c:
        resp = c.get('/api/v1/short_tandem_repeats?chrom=chr1&start=14070&end=14081&location=exact')
        data = json.loads(resp.data)
        assert len(data['results']) == 1
        assert data['results'][0]['id'] == "short_tandem_repeats_305056"
    #select_by_uids
    with app.test_client() as c:
        resp = c.get('/api/v1/short_tandem_repeats?chrom=chr1&start=14070&end=14081&location=exact&uids=305056')
        data = json.loads(resp.data)
        assert len(data['results']) == 1
        assert data['results'][0]['id'] == "short_tandem_repeats_305056"
    #select_by_sources
    with app.test_client() as c:
        resp = c.get('/api/v1/short_tandem_repeats?chrom=chr1&start=11869&end=14362&sources=gangSTR&location=overlapping')
        data = json.loads(resp.data)
        assert data['size'] == 1
    #select_by_pathogencity 
    with app.test_client() as c:
        resp = c.get('/api/v1/short_tandem_repeats/pathogenic')
        data = json.loads(resp.data)
        assert data['size'] == 26
    #select_by_motif
    with app.test_client() as c:
        resp = c.get('/api/v1/short_tandem_repeats?motif=GCT&rotation=True&chrom=chr1&start=11869&end=3004362&location=overlapping')
        data = json.loads(resp.data)
        assert data['size'] == 12
    
# def test_conservation(client): # TODO
#     #select_by_location
#     #select_by_exact_location
#     #select_by_uids
#     #select_by_sources
#     return False

# def test_repeat_mask():       # TODO
#     #select_by_location
#     #select_by_exact_location
#     #select_by_uids
#     #select_by_sources
#     return False

# # ============== gfmixin2 ============== #
# @pytest.mark.skip(reason="hg19")
def test_dna_accessibility():   
    with app.test_client() as c:        ## location exact
        resp = c.get('/api/v1/dna_accessibility?chrom=chr1&start=10410&end=10606&location=exact')
        data = json.loads(resp.data)
        assert data['size'] == 167
    with app.test_client() as c:        ## location within
        resp = c.get('/api/v1/dna_accessibility?chrom=chr1&start=10410&end=10606&location=overlapping')
        data = json.loads(resp.data)
        assert data['size'] == 168
    with app.test_client() as c:        ## uids
        resp = c.get('/api/v1/dna_accessibility?chrom=chr1&start=10410&end=10606&location=overlapping&uids=1')
        data = json.loads(resp.data)
        assert data['results'][0]['id'] == "dna_accessibility_1"
    with app.test_client() as c:        ## sources
        resp = c.get('/api/v1/dna_accessibility?chrom=chr1&start=10410&end=10606&location=overlapping&sources=ENCODE')
        data = json.loads(resp.data)
        assert data['size'] == 168
    with app.test_client() as c:        ## samples
        resp = c.get('/api/v1/dna_accessibility?chrom=chr1&start=10410&end=10606&location=overlapping&sources=ENCODE&samples=mesoderm+(heart)')
        data = json.loads(resp.data)
        assert data['size'] == 1
    with app.test_client() as c:        ## experiments 
        resp = c.get('/api/v1/dna_accessibility?chrom=chr1&start=10410&end=10606&experiments=DNase-seq&location=overlapping')
        data = json.loads(resp.data)
        assert data['size'] == 167 

# @pytest.mark.skip(reason="hg19")
def test_enhancer():            
    with app.test_client() as c:        ## location exact
        resp = c.get('/api/v1/enhancers?chrom=chr1&start=858257&end=858648&location=exact')
        data = json.loads(resp.data)
        assert data['results'][0]['id'] == "enhancers_1"
    with app.test_client() as c:        ## location within
        resp = c.get('/api/v1/enhancers?chrom=chr1&start=858257&end=858648&location=overlapping')
        data = json.loads(resp.data)
        assert data['size'] == 66
    with app.test_client() as c:        ## uids
        resp = c.get('/api/v1/enhancers?uids=1&chrom=chr1&start=858257&end=858648&location=overlapping')
        data = json.loads(resp.data)
        assert data['results'][0]['id'] == "enhancers_1"
    with app.test_client() as c:        ## sources
        resp = c.get('/api/v1/enhancers?chrom=chr1&start=858257&end=858648&sources=FANTOM5&location=overlapping')
        data = json.loads(resp.data)
        assert data['size'] == 66
    with app.test_client() as c:        ## samples
        resp = c.get('/api/v1/enhancers?chrom=chr1&start=858257&end=858648&samples=MCF-7&location=overlapping')
        data = json.loads(resp.data)
        assert data['size'] == 1
    with app.test_client() as c:        ## experiments 
        resp = c.get('/api/v1/enhancers?chrom=chr1&start=858257&end=858648&experiments=CAGE&location=overlapping')
        data = json.loads(resp.data)
        assert data['size'] == 66

# @pytest.mark.skip(reason="hg19")
def test_histone_modification():    # TODO ++
    with app.test_client() as c:        ## location exact
        resp = c.get('/api/v1/histone_modifications?chrom=chr1&start=911020&end=912066&location=exact')
        data = json.loads(resp.data)
        assert data['size'] == 37
    with app.test_client() as c:        ## location within
        resp = c.get('/api/v1/histone_modifications?chrom=chr1&start=911020&end=912066&location=overlapping')
        data = json.loads(resp.data)
        assert data['size'] == 740
    with app.test_client() as c:        ## uids
        resp = c.get('/api/v1/histone_modifications?uids=38565791,38565792&chrom=chr1&start=911020&end=912066&location=overlapping')
        data = json.loads(resp.data)
        assert data['size'] == 2
    with app.test_client() as c:        ## sources
        resp = c.get('/api/v1/histone_modifications?chrom=chr1&start=911020&end=912066&sources=ENCODE&location=overlapping')
        data = json.loads(resp.data)
        assert data['size'] == 740
    with app.test_client() as c:        ## samples
        resp = c.get('/api/v1/histone_modifications?chrom=chr1&start=911020&end=912066&samples=mucosa+(colon)&location=overlapping')
        data = json.loads(resp.data)
        assert data['size'] == 3
    with app.test_client() as c:        ## experiments 
        resp = c.get('/api/v1/histone_modifications?chrom=chr1&start=911020&end=912066&experiments=ChIP-seq&location=overlapping')
        data = json.loads(resp.data)
        assert data['size'] == 740
    with app.test_client() as c:        ## histone_type 
        resp = c.get('/api/v1/histone_modifications?chrom=chr1&start=911020&end=912066&histone_type=H3K27me3&location=overlapping')
        data = json.loads(resp.data)
        assert data['size'] == 740

# @pytest.mark.skip(reason="hg19")
def test_tad():   
    with app.test_client() as c:        ## location exact
        resp = c.get('/api/v1/tads?restriction_enzymes=HindIII&chrom=chr1&start=720001&end=3600000&location=exact')
        data = json.loads(resp.data)
        assert data['size'] == 3
    with app.test_client() as c:        ## location exact
        resp = c.get('/api/v1/tads?chrom=chr1&start=720001&end=3600000&location=exact')
        data = json.loads(resp.data)
        assert data['size'] == 3 
    with app.test_client() as c:        ## location within
        resp = c.get('/api/v1/tads?chrom=chr1&start=720001&end=3600000&location=overlapping')
        data = json.loads(resp.data)
        assert data['size'] == 71
    with app.test_client() as c:        ## uids
        resp = c.get('/api/v1/tads?uids=1&chrom=chr1&start=720001&end=3600000&location=overlapping')
        data = json.loads(resp.data)
        assert data['results'][0]['id'] == 'tads_1'
    with app.test_client() as c:        ## sources
        resp = c.get('/api/v1/tads?chrom=chr1&start=720001&end=3600000&sources=3D-genome&location=overlapping')
        data = json.loads(resp.data)
        assert data['size'] == 71
    with app.test_client() as c:        ## samples
        resp = c.get('/api/v1/tads?chrom=chr1&start=720001&end=3600000&samples=A-549&location=overlapping') ##
        data = json.loads(resp.data)
        assert data['size'] == 1
    with app.test_client() as c:        ## experiments 
        resp = c.get('/api/v1/tads?chrom=chr1&start=720001&end=3600000&experiments=Hi-C&location=overlapping') ##
        data = json.loads(resp.data)
        assert data['size'] == 36

# @pytest.mark.skip(reason="hg19")
def test_tf_binding():  
    with app.test_client() as c:        ## location exact
        resp = c.get('/api/v1/tf_binding?chrom=chr1&start=847861&end=848080&location=exact')
        data = json.loads(resp.data)
        assert data['size'] == 1
    with app.test_client() as c:        ## location within
        resp = c.get('/api/v1/tf_binding?chrom=chr1&start=847861&end=848080&location=overlapping')
        data = json.loads(resp.data)
        assert data['size'] == 117 
    with app.test_client() as c:        ## uids
        resp = c.get('/api/v1/tf_binding?uids=20852860&chrom=chr1&start=847861&end=848080&location=overlapping')
        data = json.loads(resp.data)
        assert data['results'][0]['id'] == 'tf_binding_20852860'
    with app.test_client() as c:        ## sources
        resp = c.get('/api/v1/tf_binding?chrom=chr1&start=847861&end=848080&sources=ReMap&location=overlapping')
        data = json.loads(resp.data)
        assert data['size'] == 117
    with app.test_client() as c:        ## samples
        resp = c.get('/api/v1/tf_binding?chrom=chr1&start=847861&end=848080&samples=Hep-G2&location=overlapping')
        data = json.loads(resp.data)
        assert data['size'] == 6
    with app.test_client() as c:        ## experiments 
        resp = c.get('/api/v1/tf_binding?chrom=chr1&start=847861&end=848080&experiments=ChIP-seq&location=overlapping')
        data = json.loads(resp.data)
        assert data['size'] == 117
    with app.test_client() as c:        ## tfs 
        resp = c.get('/api/v1/tf_binding?chrom=chr1&start=847861&end=848080&tfs=STAG1&location=overlapping')
        data = json.loads(resp.data)
        assert data['size'] == 2

def test_tss():         # TODO ++
    with app.test_client() as c:        ## location exact
        resp = c.get('/api/v1/tss?chrom=chr10&start=100013404&end=100013414&location=exact')
        data = json.loads(resp.data)
        assert data['results'][0]['id'] == 'transcription_start_sites_1' 
    with app.test_client() as c:        ## location within
        resp = c.get('/api/v1/tss?chrom=chr10&start=100013404&end=100053414&location=overlapping')
        data = json.loads(resp.data)
        assert data['size'] == 2 
    with app.test_client() as c:        ## uids
        resp = c.get('/api/v1/tss?chrom=chr10&start=100013404&end=100053414&location=overlapping&uids=1')
        data = json.loads(resp.data)
        assert data['size'] == 1
    with app.test_client() as c:        ## sources
        resp = c.get('/api/v1/tss?chrom=chr10&start=100013404&end=100053414&location=overlapping&sources=FANTOM5')
        data = json.loads(resp.data)
        assert data['size'] == 2 
    with app.test_client() as c:        ## experiments 
        resp = c.get('/api/v1/tss?chrom=chr10&start=100013404&end=100053414&location=overlapping&experiments=CAGE')
        data = json.loads(resp.data)
        assert data['size'] == 2 
    with app.test_client() as c:        ## genes 
        resp = c.get('/api/v1/tss?chrom=chr19&start=58858939&end=58859039&location=overlapping&genes=A1BG')
        data = json.loads(resp.data)
        assert data['size'] == 1
    with app.test_client() as c:        ## genic 
        resp = c.get('/api/v1/tss/genic')
        data = json.loads(resp.data)
        assert data['size'] == 103160

# # ============== other ============== #
# @pytest.mark.skip(reason="hg19")
def test_sample():      # TODO ++
    with app.test_client() as c:      
        resp = c.get('/api/v1/samples')
        data = json.loads(resp.data)
        assert data['size'] == 1029

# @pytest.mark.skip(reason="hg19")
def test_experiment():  # TODO ++
    with app.test_client() as c:      
        resp = c.get('/api/v1/experiments')
        data = json.loads(resp.data)
        assert data['size'] == 7

# @pytest.mark.skip(reason="hg19")
def test_source():      # TODO ++
    with app.test_client() as c:      
        resp = c.get('/api/v1/sources')
        data = json.loads(resp.data)
        assert data['size'] == 5

# @pytest.mark.skip(reason="hg19")
def test_chrom():       # TODO ++
    with app.test_client() as c:      
        resp = c.get('/api/v1/chroms')
        data = json.loads(resp.data)
        assert data['size'] == 25

# @pytest.mark.skip(reason="hg19")
def test_expression():       # TODO ++
    with app.test_client() as c:      
        resp = c.get('/api/v1/expression?samples=fibroblast+(mesenchymal+villi)&min_tpm=0&max_tpm=2')
        data = json.loads(resp.data)
        assert data['size'] == 44160
