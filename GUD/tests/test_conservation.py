import unittest
from GUD import GUDUtils
from GUD.ORM import Conservation


class ConservationTests(unittest.TestCase):
    GUDUtils.db = "test_hg38_chr22"
    db_name = GUDUtils._get_db_name()
    engine, Session = GUDUtils.get_engine_session(db_name)

    def test_select_all(self):
        session = self.Session()
        feats = Conservation.select_all(session, None).all()
        self.assertEqual(len(feats), 124988)
        session.close()
        self.engine.dispose()

    def test_select_by_overlapping_location(self):
        session = self.Session()
        feats = Conservation.select_by_location(
            session, None, "22", 29999350, 39999350, "overlapping").all()
        self.assertEqual(len(feats), 37637)
        session.close()
        self.engine.dispose()

    def test_select_by_within_location(self):
        session = self.Session()
        feats = Conservation.select_by_location(
            session, None, "22", 39978331, 39999350, "within").all()
        self.assertEqual(len(feats), 91)
        session.close()
        self.engine.dispose()

    def test_select_by_exact_location(self):
        session = self.Session()
        feats = Conservation.select_by_location(
            session, None, "22", 39999346, 39999350, "exact").all()
        self.assertEqual(len(feats), 1)
        session.close()
        self.engine.dispose()

    def test_select_by_uids(self):
        session = self.Session()
        feats = Conservation.select_by_uids(
            session, None, [3431560, 3431571]).all()
        self.assertEqual(len(feats), 2)
        session.close()
        self.engine.dispose()

    def test_select_by_sources(self):
        session = self.Session()
        feats = Conservation.select_by_sources(session, None , ["phastConsElements100way"]).all()
        self.assertEqual(len(feats), 124988)
        session.close()
        self.engine.dispose()

    def test_compound_select(self):
        session = self.Session()
        feats = Conservation.select_by_location(
            session, None, "22",  39978331, 39999350, "within")
        feats = Conservation.select_by_uids(session, feats, [4434805, 4434772]).all()
        self.assertEqual(len(feats), 2)
        session.close()
        self.engine.dispose()

    def test_is_unique(self):
        session = self.Session()
        feats = Conservation.is_unique(session, 3474989,2)
        self.assertFalse(feats)
        session.close()
        self.engine.dispose()

    def test_as_genomic_feature(self):
        session = self.Session()
        feats = Conservation.select_by_uids(
            session, None, [4434805]).first()
        gf = Conservation.as_genomic_feature(feats)
        self.assertEqual(gf.chrom, "22")
        self.assertEqual(gf.start, 39999346)
        self.assertEqual(gf.end, 39999350)
        self.assertEqual(gf.score, 1)
        self.assertEqual(gf.strand, 0)
        self.assertEqual(gf.id, "conservation_4434805")
        self.assertEqual(gf.qualifiers["source"], "phastConsElements100way")
        self.assertEqual(gf.qualifiers["uid"], 4434805)
        session.close()
        self.engine.dispose()


if __name__ == '__main__':
    unittest.main()
