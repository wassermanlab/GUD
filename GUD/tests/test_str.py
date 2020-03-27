import unittest
from GUD import GUDUtils
from GUD.ORM import ShortTandemRepeat


class ShortTandemRepeatTests(unittest.TestCase):
    GUDUtils.db = "test_hg38_chr22"
    db_name = GUDUtils._get_db_name()
    engine, Session = GUDUtils.get_engine_session(db_name)

    def test_select_all(self):
        session = self.Session()
        feats = ShortTandemRepeat.select_all(session, None).all()
        self.assertEqual(len(feats), 11904)
        session.close()
        self.engine.dispose()

    def test_select_by_overlapping_location(self):
        session = self.Session()
        feats = ShortTandemRepeat.select_by_location(
            session, None, "22", 50008263, 50808291, "overlapping").all()
        self.assertEqual(len(feats), 262)
        session.close()
        self.engine.dispose()

    def test_select_by_within_location(self):
        session = self.Session()
        feats = ShortTandemRepeat.select_by_location(
            session, None, "22", 50008263, 50808291, "within").all()
        self.assertEqual(len(feats), 262)
        session.close()
        self.engine.dispose()

    def test_select_by_exact_location(self):
        session = self.Session()
        feats = ShortTandemRepeat.select_by_location(
            session, None, "22", 50808263, 50808291, "exact").all()
        self.assertEqual(len(feats), 1)
        session.close()
        self.engine.dispose()

    def test_select_by_uids(self):
        session = self.Session()
        feats = ShortTandemRepeat.select_by_uids(
            session, None, [318236, 318244]).all()
        self.assertEqual(len(feats), 2)
        session.close()
        self.engine.dispose()

    def test_select_by_sources(self):
        session = self.Session()
        feats = ShortTandemRepeat.select_by_sources(session, None , ["GangSTR"]).all()
        self.assertEqual(len(feats), 11903)
        session.close()
        self.engine.dispose()

    def test_compound_select(self):
        session = self.Session()
        feats = ShortTandemRepeat.select_by_location(
            session, None, "22", 50008263, 50808291, "within")
        feats = ShortTandemRepeat.select_by_uids(session, feats, [318236, 318244]).all()
        self.assertEqual(len(feats), 2)
        session.close()
        self.engine.dispose()

    def test_is_unique(self):
        session = self.Session()
        feats = ShortTandemRepeat.is_unique(session, 4556568, 4, 266)
        self.assertFalse(feats)
        session.close()
        self.engine.dispose()

    def test_as_genomic_feature(self):
        session = self.Session()
        feats = ShortTandemRepeat.select_by_uids(
            session, None, [15]).first()
        gf = ShortTandemRepeat.as_genomic_feature(feats)
        self.assertEqual(gf.chrom, "22")
        self.assertEqual(gf.start, 45795354)
        self.assertEqual(gf.end, 45795424)
        self.assertEqual(gf.strand, 0)
        self.assertEqual(gf.id, "short_tandem_repeats_15")
        self.assertEqual(gf.qualifiers["motif"], "ATTCT")
        self.assertEqual(gf.qualifiers["pathogenicity"], 266)
        self.assertEqual(gf.qualifiers["source"], "Richmond_github_20191007")
        self.assertEqual(gf.qualifiers["uid"], 15)
        session.close()
        self.engine.dispose()

    def test_select_by_pathogenicity(self):
        session = self.Session()
        feats = ShortTandemRepeat.select_by_pathogenicity(session, None ).all()
        self.assertEqual(len(feats), 1)
        session.close()
        self.engine.dispose()

    def test_select_by_motif(self):
        session = self.Session()
        feats = ShortTandemRepeat.select_by_motif(session, "AAT", None, False).all()
        feats_rotations = ShortTandemRepeat.select_by_motif(session, "AAT", None, True).all()
        self.assertEqual(len(feats), 205)
        self.assertEqual(len(feats_rotations), 319)
        session.close()
        self.engine.dispose()


if __name__ == '__main__':
    unittest.main()
