import unittest
import os, sys
from GUD.api import app
from flask import json
from time import sleep

class STRTests(unittest.TestCase):
    def setUp(self):
        # creates a test client
        self.app = app.test_client()
        # propagate the exceptions to the test client
        self.app.testing = True 

    def tearDown(self):
        sleep(1)
        pass 
    
    def test_select_all(self):
        resp = self.app.get('/api/v1/test_hg38_chr22/short_tandem_repeats')
        data = json.loads(resp.data)
        self.assertEqual(data["size"], 11904)

    def test_select_by_overlapping_location(self):
        resp = self.app.get('/api/v1/test_hg38_chr22/short_tandem_repeats?chrom=22&start=50008264&end=50808291&location=overlapping')
        data = json.loads(resp.data)
        self.assertEqual(data["size"], 262)

    def test_select_by_within_location(self):
        resp = self.app.get('/api/v1/test_hg38_chr22/short_tandem_repeats?chrom=22&start=50008264&end=50808291&location=within')
        data = json.loads(resp.data)
        self.assertEqual(data["size"], 262)

    def test_select_by_exact_location(self):
        resp = self.app.get('/api/v1/test_hg38_chr22/short_tandem_repeats?chrom=22&start=50808264&end=50808291&location=exact')
        data = json.loads(resp.data)
        self.assertEqual(data["size"], 1)

    def test_select_by_uids(self):
        resp = self.app.get('/api/v1/test_hg38_chr22/short_tandem_repeats?uids=318236,318244')
        data = json.loads(resp.data)
        self.assertEqual(data["size"], 2)

    def test_select_by_sources(self):
        resp = self.app.get('/api/v1/test_hg38_chr22/short_tandem_repeats?sources=GangSTR')
        data = json.loads(resp.data)
        self.assertEqual(data["size"], 11903)

    def test_compound_select(self):
        resp = self.app.get('/api/v1/test_hg38_chr22/short_tandem_repeats?chrom=22&start=50008263&end=50808291&location=within&uids=318236,318244')
        data = json.loads(resp.data)
        self.assertEqual(data["size"], 2)

    def test_select_by_pathogenicity(self):
        resp = self.app.get('/api/v1/test_hg38_chr22/short_tandem_repeats/pathogenic')
        data = json.loads(resp.data)
        self.assertEqual(data["size"], 1)

    def test_select_by_motif(self):
        resp = self.app.get('/api/v1/test_hg38_chr22/short_tandem_repeats?motif=AAT')
        data = json.loads(resp.data)
        self.assertEqual(data["size"], 205)

    def test_select_by_motif_rotations(self):
        resp = self.app.get('/api/v1/test_hg38_chr22/short_tandem_repeats?motif=AAT&rotation=True')
        data = json.loads(resp.data)
        self.assertEqual(data["size"], 319)

    def test_formatting(self):
        resp = self.app.get('/api/v1/test_hg38_chr22/short_tandem_repeats?chrom=22&start=50008263&end=50808291&location=within&uids=318236,318244')
        data = json.loads(resp.data)
        self.assertEqual(data["results"][0]["chrom"], "22")
        self.assertEqual(data["results"][0]["start"], 50797605)
        self.assertEqual(data["results"][0]["end"], 50797625)
        self.assertEqual(data["results"][0]["id"], "short_tandem_repeats_318236")

    def test_pagination(self):
        resp = self.app.get('/api/v1/test_hg38_chr22/short_tandem_repeats')
        data = json.loads(resp.data)
        self.assertEqual(data["size"], 11904)
        self.assertEqual(len(data["results"]), 20)
        self.assertEqual(data["next"], 'http://localhost/api/v1/test_hg38_chr22/short_tandem_repeats?page=2')


if __name__ == '__main__':
    unittest.main()
