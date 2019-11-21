import unittest
import os, sys
from GUD.api import app
from flask import json
from time import sleep

class CNVTests(unittest.TestCase):
    def setUp(self):
        # creates a test client
        self.app = app.test_client()
        # propagate the exceptions to the test client
        self.app.testing = True 

    def tearDown(self):
        sleep(1)
        pass 

    def test_select_all(self):
        resp = self.app.get('/api/v1/test_hg38_chr22/copy_number_variants')
        data = json.loads(resp.data)
        self.assertEqual(data["next"], "http://localhost/api/v1/test_hg38_chr22/copy_number_variants?last_uid=2788")

    def test_select_by_overlapping_location(self):
        resp = self.app.get('/api/v1/test_hg38_chr22/copy_number_variants?chrom=22&start=49199842&end=49787792&location=overlapping')
        data = json.loads(resp.data)
        self.assertEqual(data["next"], "http://localhost/api/v1/test_hg38_chr22/copy_number_variants?chrom=22&start=49199842&end=49787792&location=overlapping&last_uid=6068")

    def test_select_by_within_location(self):
        resp = self.app.get('/api/v1/test_hg38_chr22/copy_number_variants?chrom=22&start=49199842&end=49787792&location=within')
        data = json.loads(resp.data)
        self.assertEqual(len(data["results"]), 1)

    def test_select_by_exact_location(self):
        resp = self.app.get('/api/v1/test_hg38_chr22/copy_number_variants?chrom=22&start=49199842&end=49787792&location=exact')
        data = json.loads(resp.data)
        self.assertEqual(len(data["results"]), 1)

    def test_select_by_uids(self):
        resp = self.app.get('/api/v1/test_hg38_chr22/copy_number_variants?uids=11849,11843')
        data = json.loads(resp.data)
        self.assertEqual(len(data["results"]), 2)

    def test_select_by_sources(self):
        resp = self.app.get('/api/v1/test_hg38_chr22/copy_number_variants?sources=dbVar')
        data = json.loads(resp.data)
        self.assertEqual(data["next"], "http://localhost/api/v1/test_hg38_chr22/copy_number_variants?sources=dbVar&last_uid=2788")

    def test_compound_select(self):
        resp = self.app.get('/api/v1/test_hg38_chr22/copy_number_variants?chrom=22&start=49199841&end=49787792&location=overlapping&uids=11849,11837')
        data = json.loads(resp.data)
        self.assertEqual(len(data["results"]), 2)

    def select_by_clinical_assertion(self):
        resp = self.app.get('/api/v1/test_hg38_chr22/copy_number_variants?clinical_assertion=Pathogenic')
        data = json.loads(resp.data)
        self.assertEqual(data["next"], "http://localhost/api/v1/test_hg38_chr22/copy_number_variants?clinical_assertion=Pathogenic&last_uid=2788")

    def test_select_by_clinvar_accession(self):
        resp = self.app.get('/api/v1/test_hg38_chr22/copy_number_variants?clinvar_accession=SCV000080129')
        data = json.loads(resp.data)
        self.assertEqual(len(data["results"]), 1)
    
    def test_selec_by_dnvar_accession(self):
        resp = self.app.get('/api/v1/test_hg38_chr22/copy_number_variants?dbVar_accession=nssv15120684')
        data = json.loads(resp.data)
        self.assertEqual(len(data["results"]), 1)

if __name__ == '__main__':
    unittest.main()
