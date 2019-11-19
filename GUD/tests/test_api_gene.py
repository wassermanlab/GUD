import unittest
import os, sys
from GUD.api import app
from flask import json
from time import sleep

class GeneTests(unittest.TestCase):
    def setUp(self):
        # creates a test client
        self.app = app.test_client()
        # propagate the exceptions to the test client
        self.app.testing = True 

    def tearDown(self):
        sleep(1)
        pass 

    def test_select_all(self):
        resp = self.app.get('/api/v1/test_hg38_chr22/genes')
        data = json.loads(resp.data)
        self.assertEqual(data["size"], 1554)

    def test_select_by_overlapping_location(self):
        resp = self.app.get('/api/v1/test_hg38_chr22/genes?chrom=22&start=1&end=15000000&location=overlapping')
        data = json.loads(resp.data)
        self.assertEqual(data["size"], 2)

    def test_select_by_within_location(self):
        resp = self.app.get('/api/v1/test_hg38_chr22/genes?chrom=22&start=1&end=15000000&location=within')
        data = json.loads(resp.data)
        self.assertEqual(data["size"], 2)

    def test_select_by_exact_location(self):
        resp = self.app.get('/api/v1/test_hg38_chr22/genes?chrom=22&start=10940597&end=10961529&location=exact')
        data = json.loads(resp.data)
        self.assertEqual(data["size"], 1)

    def test_select_by_uids(self):
        resp = self.app.get('/api/v1/test_hg38_chr22/genes?uids=34740')
        data = json.loads(resp.data)
        self.assertEqual(data["size"], 1)

    def test_select_by_sources(self):
        resp = self.app.get('/api/v1/test_hg38_chr22/genes?sources=refGene')
        data = json.loads(resp.data)
        self.assertEqual(data["size"], 1554)

    def test_compound_select(self):
        resp = self.app.get('/api/v1/test_hg38_chr22/genes?chrom=22&start=1&end=150000001&location=within&uids=34740')
        data = json.loads(resp.data)
        self.assertEqual(data["size"], 1)

    def test_select_by_names(self):
        resp = self.app.get('/api/v1/test_hg38_chr22/genes?names=YDJC,XBP1')
        data = json.loads(resp.data)
        self.assertEqual(data["size"], 8)
    
    def test_get_all_gene_symbols(self):
        resp = self.app.get('/api/v1/test_hg38_chr22/genes/symbols')
        data = json.loads(resp.data)
        self.assertEqual(data["size"], 644)

if __name__ == '__main__':
    unittest.main()
