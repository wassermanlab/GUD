import unittest
import os, sys
from GUD.api import app
import json

class ChromTests(unittest.TestCase):
    def setUp(self):
        # creates a test client
        self.app = app.test_client()
        # propagate the exceptions to the test client
        self.app.testing = True 

    def tearDown(self):
        pass 

    def test_all_chroms(self):
        pass

    def test_select_chroms(self):
        pass
    
    def test_select_sizes(self):
        pass

if __name__ == '__main__':
    unittest.main()
