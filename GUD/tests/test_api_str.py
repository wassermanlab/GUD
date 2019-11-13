import unittest
import os, sys
from GUD.api import app
import json

class STRTests(unittest.TestCase):
    def setUp(self):
        # creates a test client
        self.app = app.test_client()
        # propagate the exceptions to the test client
        self.app.testing = True 

    def tearDown(self):
        pass 

    def test_select_all(self):
        pass

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

    def test_compound_select(self):
        pass

    def test_select_by_pathogenicity(self):
        pass

    def test_select_by_motif(self):
        pass

    def test_select_by_motif_rotations(self):
        pass


if __name__ == '__main__':
    unittest.main()
