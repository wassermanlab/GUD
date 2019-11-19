import unittest
import os, sys
from GUD.api import app
from flask import json
from time import sleep


class RegionTests(unittest.TestCase):
    def setUp(self):
        # creates a test client
        self.app = app.test_client()
        # propagate the exceptions to the test client
        self.app.testing = True 

    def tearDown(self):
        sleep(1)
        pass

    def test_select_all(self):
        resp = self.app.get('/api/v1/test_hg38_chr22/clinvar')
        data = json.loads(resp.data)
        self.assertEqual(data["size"], 7943)

    def test_select_by_overlapping_location(self):
        resp = self.app.get('/api/v1/test_hg38_chr22/clinvar?chrom=22&start=50644827&end=50744826&location=overlapping')
        data = json.loads(resp.data)
        self.assertEqual(data["size"], 156)

    def test_select_by_within_location(self):
        resp = self.app.get('/api/v1/test_hg38_chr22/clinvar?chrom=22&start=50644827&end=50744826&location=within')
        data = json.loads(resp.data)
        self.assertEqual(data["size"], 156)

    def test_select_by_exact_location(self):
        resp = self.app.get('/api/v1/test_hg38_chr22/clinvar?chrom=22&start=50744827&end=50744827&location=exact')
        data = json.loads(resp.data)
        self.assertEqual(data["size"], 1)

    def test_select_by_uids(self):
        resp = self.app.get('/api/v1/test_hg38_chr22/clinvar?uids=246848,246854')
        data = json.loads(resp.data)
        self.assertEqual(data["size"], 2)

    def test_select_by_sources(self):
        resp = self.app.get('/api/v1/test_hg38_chr22/clinvar?sources=ClinVar')
        data = json.loads(resp.data)
        self.assertEqual(data["size"], 7943)

    def test_compound_select(self):
        resp = self.app.get('/api/v1/test_hg38_chr22/clinvar?chrom=22&start=50644826&end=50744826&location=within&uids=246848,246854')
        data = json.loads(resp.data)
        self.assertEqual(data["size"], 2)

    def test_select_by_clinvarID(self):
        resp = self.app.get('/api/v1/test_hg38_chr22/clinvar?clinvar_ids=51')
        data = json.loads(resp.data)
        self.assertEqual(data["size"], 1)

if __name__ == '__main__':
    unittest.main()
