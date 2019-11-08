import unittest
from GUD import GUDUtils
from GUD.ORM import Chrom

class RegionTests(unittest.TestCase):
    GUDUtils.db = "test_hg38_chr22"
    db_name = GUDUtils._get_db_name()
    engine, Session = GUDUtils.get_engine_session(db_name)

    def test_all_chroms(self):
        session = self.Session()
        feats = Chrom.select_all_chroms(session).all()
        self.assertEqual(len(feats), 25)
        session.close()
        self.engine.dispose() 

    def test_select_chroms(self):
        session = self.Session()
        feats = Chrom.select_by_chroms(session, ["1", "5"])
        self.assertEqual(len(feats), 2)
        session.close()
        self.engine.dispose() 
    
    def test_select_sizes(self):
        session = self.Session()
        feats = Chrom.chrom_sizes(session, ["1", "5"])
        self.assertEqual(len(feats), 2)
        session.close()
        self.engine.dispose() 

if __name__ == '__main__':
    unittest.main()
