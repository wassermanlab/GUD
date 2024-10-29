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
        self.assertEqual(data["next"], "http://localhost/api/v1/test_hg38_chr22/short_tandem_repeats?last_uid=258787")

    def test_select_by_overlapping_location(self):
        resp = self.app.get('/api/v1/test_hg38_chr22/short_tandem_repeats?chrom=22&start=50008264&end=50808291&location=overlapping')
        data = json.loads(resp.data)
        self.assertEqual(data["next"], "http://localhost/api/v1/test_hg38_chr22/short_tandem_repeats?chrom=22&start=50008264&end=50808291&location=overlapping&last_uid=317019")

    def test_select_by_within_location(self):
        resp = self.app.get('/api/v1/test_hg38_chr22/short_tandem_repeats?chrom=22&start=50008264&end=50808291&location=within')
        data = json.loads(resp.data)
        self.assertEqual(data["next"], "http://localhost/api/v1/test_hg38_chr22/short_tandem_repeats?chrom=22&start=50008264&end=50808291&location=within&last_uid=317019")

    def test_select_by_exact_location(self):
        resp = self.app.get('/api/v1/test_hg38_chr22/short_tandem_repeats?chrom=22&start=50808264&end=50808291&location=exact')
        data = json.loads(resp.data)
        self.assertEqual(len(data["results"]), 1)

    def test_select_by_uids(self):
        resp = self.app.get('/api/v1/test_hg38_chr22/short_tandem_repeats?uids=318236,318244')
        data = json.loads(resp.data)
        self.assertEqual(len(data["results"]), 2)

    def test_select_by_sources(self):
        resp = self.app.get('/api/v1/test_hg38_chr22/short_tandem_repeats?sources=GangSTR')
        data = json.loads(resp.data)
        self.assertEqual(data["next"], "http://localhost/api/v1/test_hg38_chr22/short_tandem_repeats?sources=GangSTR&last_uid=258791")

    def test_compound_select(self):
        resp = self.app.get('/api/v1/test_hg38_chr22/short_tandem_repeats?chrom=22&start=50008263&end=50808291&location=within&uids=318236,318244')
        data = json.loads(resp.data)
        self.assertEqual(len(data["results"]), 2)

    def test_select_by_motif(self):
        resp = self.app.get('/api/v1/test_hg38_chr22/short_tandem_repeats?motif=AAT')
        data = json.loads(resp.data)
        self.assertEqual(data["next"], "http://localhost/api/v1/test_hg38_chr22/short_tandem_repeats?motif=AAT&last_uid=265235")

    def test_select_by_motif_rotations(self):
        resp = self.app.get('/api/v1/test_hg38_chr22/short_tandem_repeats?motif=AAT&rotation=True')
        data = json.loads(resp.data)
        self.assertEqual(data["next"], "http://localhost/api/v1/test_hg38_chr22/short_tandem_repeats?motif=AAT&rotation=True&last_uid=263192")

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
        self.assertEqual(len(data["results"]), 20)
        self.assertEqual(data["next"], 'http://localhost/api/v1/test_hg38_chr22/short_tandem_repeats?last_uid=258787')


if __name__ == '__main__':
    unittest.main()
