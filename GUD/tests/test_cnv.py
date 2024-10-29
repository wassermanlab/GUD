import unittest
from GUD import GUDUtils
from GUD.ORM import CNV


class CNVTests(unittest.TestCase):
    GUDUtils.db = "test_hg38_chr22"
    db_name = GUDUtils._get_db_name()
    engine, Session = GUDUtils.get_engine_session(db_name)

    def test_select_all(self):
        session = self.Session()
        feats = CNV.select_all(session, None).all()
        self.assertEqual(len(feats), 603)
        session.close()
        self.engine.dispose()

    def test_select_by_overlapping_location(self):
        session = self.Session()
        feats = CNV.select_by_location(
            session, None, "22", 49199841, 49787792, "overlapping").all()
        self.assertEqual(len(feats), 89)
        session.close()
        self.engine.dispose()

    def test_select_by_within_location(self):
        session = self.Session()
        feats = CNV.select_by_location(
            session, None, "22", 49199841, 49787792, "within").all()
        self.assertEqual(len(feats), 1)
        session.close()
        self.engine.dispose()

    def test_select_by_exact_location(self):
        session = self.Session()
        feats = CNV.select_by_location(
            session, None, "22", 49199841, 49787792, "exact").all()
        self.assertEqual(len(feats), 1)
        session.close()
        self.engine.dispose()

    def test_select_by_uids(self):
        session = self.Session()
        feats = CNV.select_by_uids(
            session, None, [11849, 11843]).all()
        self.assertEqual(len(feats), 2)
        session.close()
        self.engine.dispose()

    def test_select_by_sources(self):
        session = self.Session()
        feats = CNV.select_by_sources(session, None , ["dbVar"]).all()
        self.assertEqual(len(feats), 603)
        session.close()
        self.engine.dispose()

    def test_compound_select(self):
        session = self.Session()
        feats = CNV.select_by_location(
            session, None, "22", 49199841, 49787792, "overlapping")
        feats = CNV.select_by_uids(session, feats, [11849, 11837]).all()
        self.assertEqual(len(feats), 2)
        session.close()
        self.engine.dispose()

    def test_is_unique(self):
        session = self.Session()
        feats = CNV.is_unique(session, 4367620, 3, -1)
        self.assertFalse(feats)
        session.close()
        self.engine.dispose()

    def test_as_genomic_feature(self):
        session = self.Session()
        feats = CNV.select_by_uids(
            session, None, [2568]).first()
        gf = CNV.as_genomic_feature(feats)
        self.assertEqual(gf.chrom, "22")
        self.assertEqual(gf.start, 10761695)
        self.assertEqual(gf.end, 16243442)
        self.assertEqual(gf.strand, 0)
        self.assertEqual(gf.id, "copy_number_variants_2568")
        self.assertEqual(gf.qualifiers["clinical_assertion"], b"Pathogenic")
        self.assertEqual(gf.qualifiers["clinvar_accession"], b"RCV000052775.5,SCV000080129")
        self.assertEqual(gf.qualifiers["dbVar_accession"], b"nssv15120684,nssv577812,essv6983991,essv6989290")
        self.assertEqual(gf.qualifiers["source"], "dbVar")
        self.assertEqual(gf.qualifiers["uid"], 2568)
        session.close()
        self.engine.dispose()

    def select_by_clinical_assertion(self):
        session = self.Session()
        feats = CNV.select_by_sources(session, None , "Pathogenic").all()
        self.assertEqual(len(feats), 581)
        session.close()
        self.engine.dispose()

    def test_select_by_clinvar_accession(self):
        session = self.Session()
        feats = CNV.select_by_clinvar_accession(session, None , "SCV000080129").all()
        self.assertEqual(len(feats), 1)
        session.close()
        self.engine.dispose()
    
    def test_selec_by_dnvar_accession(self):
        session = self.Session()
        feats = CNV.select_by_dbvar_accession(session, None , "nssv15120684").all()
        self.assertEqual(len(feats), 1)
        session.close()
        self.engine.dispose()
    

if __name__ == '__main__':
    unittest.main()
