import unittest
import os, sys
from GUD.api import app
from flask import json

class GeneTests(unittest.TestCase):
    def setUp(self):
        # creates a test client
        self.app = app.test_client()
        # propagate the exceptions to the test client
        self.app.testing = True 

    def tearDown(self):
        pass 

    def test_select_all(self):
        resp = self.app.get('/api/v1/test_hg38_chr22/chroms')
        data = json.loads(resp.data)
        self.assertEqual(len(data["results"]), 20)
        self.assertEqual(data["size"], 25)
        self.assertEqual(data["results"][0]["chrom"], "1")

    def test_select_by_overlapping_location(self):
        pass

    def test_select_by_within_location(self):
        pass

    def test_select_by_exact_location(self):
        pass

    def test_select_by_uids(self):
        pass

    def test_select_by_sources(self):
        pass

    def test_select_by_names(self):
        pass
    
    def test_get_all_gene_symbols(self):
        pass
    
    def test_compound_select(self):
        pass

if __name__ == '__main__':
    unittest.main()
