import unittest
from GUD import GUDUtils
from GUD.ORM import Gene


class RegionTests(unittest.TestCase):
    GUDUtils.db = "test_hg38_chr22"
    db_name = GUDUtils._get_db_name()
    engine, Session = GUDUtils.get_engine_session(db_name)

    def test_select_all(self):
        session = self.Session()
        feats = Gene.select_all(session, None).all()
        self.assertEqual(len(feats), 1554)
        session.close()
        self.engine.dispose()

    def test_select_by_overlapping_location(self):
        session = self.Session()
        feats = Gene.select_by_location(
            session, None, "22", 1, 15000000, "overlapping").all()
        self.assertEqual(len(feats), 2)
        session.close()
        self.engine.dispose()

    def test_select_by_within_location(self):
        session = self.Session()
        feats = Gene.select_by_location(
            session, None, "22", 1, 15000000, "within").all()
        self.assertEqual(len(feats), 2)
        session.close()
        self.engine.dispose()

    def test_select_by_exact_location(self):
        session = self.Session()
        feats = Gene.select_by_location(
            session, None, "22", 10940596, 10961529, "exact").all()
        self.assertEqual(len(feats), 1)
        session.close()
        self.engine.dispose()

    def test_select_by_uids(self):
        session = self.Session()
        feats = Gene.select_by_uids(
            session, None, [34740]).all()
        self.assertEqual(len(feats), 1)
        session.close()
        self.engine.dispose()

    def test_select_by_sources(self):
        session = self.Session()
        feats = Gene.select_by_sources(session, None , ["refGene"]).all()
        self.assertEqual(len(feats), 1554)
        session.close()
        self.engine.dispose()

    def test_compound_select(self):
        session = self.Session()
        feats = Gene.select_by_location(
            session, None, "22", 1, 15000000, "within")
        feats = Gene.select_by_uids(session, feats, [34740]).all()
        self.assertEqual(len(feats), 1)
        session.close()
        self.engine.dispose()

    def test_is_unique(self):
        session = self.Session()
        feats = Gene.is_unique(session, 17658, 1, "NM_013416")
        self.assertFalse(feats)
        session.close()
        self.engine.dispose()

    def test_as_genomic_feature(self):
        session = self.Session()
        feats = Gene.select_by_uids(
            session, None, [26950]).first()
        gf = Gene.as_genomic_feature(feats)
        self.assertEqual(gf.chrom, "22")
        self.assertEqual(gf.start, 36860987)
        self.assertEqual(gf.end, 36878017)
        self.assertEqual(gf.strand, 1)
        self.assertEqual(gf.id, "genes_26950")
        self.assertEqual(gf.qualifiers["accession_number"], "NM_013416")
        self.assertEqual(gf.qualifiers["gene_symbol"], "NCF4")
        self.assertEqual(gf.qualifiers["cdsStart"], 36861171)
        self.assertEqual(gf.qualifiers["cdsEnd"], 36876072)
        self.assertEqual(gf.qualifiers["exonStarts"], b"36860987,36864044,36864918,36867391,36870414,36871651,36872326,36875652,36877627,")
        self.assertEqual(gf.qualifiers["exonEnds"], b"36861203,36864129,36865072,36867462,36870542,36871709,36872425,36876094,36878017,")
        self.assertEqual(gf.qualifiers["source"], "refGene")
        self.assertEqual(gf.qualifiers["uid"], 26950)
        session.close()
        self.engine.dispose()

    def test_select_by_names(self):
        session = self.Session()
        feats = Gene.select_by_names(session, None ,["YDJC", "XBP1"]).all()
        self.assertEqual(len(feats), 8)
        session.close()
        self.engine.dispose()
    
    def test_get_all_gene_symbols(self):
        session = self.Session()
        feats = Gene.get_all_gene_symbols(session).all()
        self.assertEqual(len(feats), 644)
        session.close()
        self.engine.dispose()

if __name__ == '__main__':
    unittest.main()
