import unittest
import os, sys
from GUD.api import app
from flask import json
from time import sleep

class ConservationTests(unittest.TestCase):
    def setUp(self):
        # creates a test client
        self.app = app.test_client()
        # propagate the exceptions to the test client
        self.app.testing = True 

    def tearDown(self):
        sleep(1)
        pass 

    def test_select_all(self):
        resp = self.app.get('/api/v1/test_hg38_chr22/conservation')
        data = json.loads(resp.data)
        self.assertEqual(data["size"], 124988)

    def test_select_by_overlapping_location(self):
        resp = self.app.get('/api/v1/test_hg38_chr22/conservation?chrom=22&start=29999351&end=39999350&location=overlapping')
        data = json.loads(resp.data)
        self.assertEqual(data["size"], 37637)

    def test_select_by_within_location(self):
        resp = self.app.get('/api/v1/test_hg38_chr22/conservation?chrom=22&start=39978332&end=39999350&location=within')
        data = json.loads(resp.data)
        self.assertEqual(data["size"], 91)

    def test_select_by_exact_location(self):
        resp = self.app.get('/api/v1/test_hg38_chr22/conservation?chrom=22&start=39999347&end=39999350&location=exact')
        data = json.loads(resp.data)
        self.assertEqual(data["size"], 1)

    def test_select_by_uids(self):
        resp = self.app.get('/api/v1/test_hg38_chr22/conservation?uids=3431560,3431571')
        data = json.loads(resp.data)
        self.assertEqual(data["size"], 2)

    def test_select_by_sources(self):
        resp = self.app.get('/api/v1/test_hg38_chr22/conservation?sources=phastConsElements100way')
        data = json.loads(resp.data)
        self.assertEqual(data["size"], 124988)

    def test_compound_select(self):
        resp = self.app.get('/api/v1/test_hg38_chr22/conservation?chrom=22&start=39978332&end=39999350&location=within&uids=4434805,4434772')
        data = json.loads(resp.data)
        self.assertEqual(data["size"], 2)


if __name__ == '__main__':
    unittest.main()
