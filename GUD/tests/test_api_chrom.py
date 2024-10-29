import unittest
import os, sys
from GUD.api import app
from flask import json

class ChromTests(unittest.TestCase):
    def setUp(self):
        # creates a test client
        self.app = app.test_client()
        # propagate the exceptions to the test client
        self.app.testing = True 

    def tearDown(self):
        pass 

    def test_all_chroms(self):
        resp = self.app.get('/api/v1/grch37/chroms')
        data = json.loads(resp.data)
        print(data)
        # self.assertEqual(len(data["results"]), 25)
        # self.assertEqual(data["results"][0]["chrom"], "1")

if __name__ == '__main__':
    unittest.main()
