import unittest
from GUD import GUDUtils
from GUD.ORM import ClinVar


class RegionTests(unittest.TestCase):
    GUDUtils.db = "test_hg38_chr22"
    db_name = GUDUtils._get_db_name()
    engine, Session = GUDUtils.get_engine_session(db_name)

    def test_select_all(self):
        session = self.Session()
        feats = ClinVar.select_all(session, None).all()
        self.assertEqual(len(feats), 7943)
        session.close()
        self.engine.dispose()

    def test_select_by_overlapping_location(self):
        session = self.Session()
        feats = ClinVar.select_by_location(
            session, None, "22", 50644826, 50744826, "overlapping").all()
        self.assertEqual(len(feats), 156)
        session.close()
        self.engine.dispose()

    def test_select_by_within_location(self):
        session = self.Session()
        feats = ClinVar.select_by_location(
            session, None, "22", 50644826, 50744826, "within").all()
        self.assertEqual(len(feats), 156)
        session.close()
        self.engine.dispose()

    def test_select_by_exact_location(self):
        session = self.Session()
        feats = ClinVar.select_by_location(
            session, None, "22", 50744826, 50744827, "exact").all()
        self.assertEqual(len(feats), 1)
        session.close()
        self.engine.dispose()

    def test_select_by_uids(self):
        session = self.Session()
        feats = ClinVar.select_by_uids(
            session, None, [246848, 246854]).all()
        self.assertEqual(len(feats), 2)
        session.close()
        self.engine.dispose()

    def test_select_by_sources(self):
        session = self.Session()
        feats = ClinVar.select_by_sources(session, None , ["ClinVar"]).all()
        self.assertEqual(len(feats), 7943)
        session.close()
        self.engine.dispose()

    def test_select_by_clinvarID(self):
        session = self.Session()
        feats = ClinVar.select_by_clinvarID(session, None ,[51]).all()
        self.assertEqual(len(feats), 1)
        session.close()
        self.engine.dispose()

    def test_compound_select(self):
        session = self.Session()
        feats = ClinVar.select_by_location(
            session, None, "22", 50644826, 50744826, "within")
        feats = ClinVar.select_by_uids(session, feats, [246854, 246848]).all()
        self.assertEqual(len(feats), 2)
        session.close()
        self.engine.dispose()

    def test_is_unique(self):
        session = self.Session()
        feats = ClinVar.is_unique(session, 51)
        self.assertFalse(feats)
        session.close()
        self.engine.dispose()

    def test_as_genomic_feature(self):
        session = self.Session()
        feats = ClinVar.select_by_uids(
            session, None, [246854]).first()
        gf = ClinVar.as_genomic_feature(feats)
        self.assertEqual(gf.chrom, "22")
        self.assertEqual(gf.start, 50744056)
        self.assertEqual(gf.end, 50744057)
        self.assertEqual(gf.strand, 0)
        self.assertEqual(gf.id, "clinvar_246854")
        self.assertEqual(gf.qualifiers["ANN_Annotation"], "splice_region_variant&intron_variant")
        self.assertEqual(gf.qualifiers["ANN_Annotation_Impact"], "LOW")
        self.assertEqual(gf.qualifiers["ANN_Feature_ID"], "ENST00000216139.9")
        self.assertEqual(gf.qualifiers["ANN_Feature_Type"], "transcript")
        self.assertEqual(gf.qualifiers["ANN_Gene_ID"], "ENSG00000100312")
        self.assertEqual(gf.qualifiers["ANN_Gene_Name"], "ACR")
        self.assertEqual(gf.qualifiers["CADD"], 6.913)
        self.assertEqual(gf.qualifiers["CLNDISDB"], "MedGen:CN169374")
        self.assertEqual(gf.qualifiers["CLNDN"], "not_specified")
        self.assertEqual(gf.qualifiers["CLNSIG"], "Benign")
        self.assertEqual(gf.qualifiers["alt"], "G")
        self.assertEqual(gf.qualifiers["clinvarID"], 402335)
        self.assertEqual(gf.qualifiers["gnomad_exome_af_global"], None)
        self.assertEqual(gf.qualifiers["gnomad_exome_hom_global"], None)
        self.assertEqual(gf.qualifiers["gnomad_genome_af_global"], None)
        self.assertEqual(gf.qualifiers["gnomad_genome_hom_global"], None)
        self.assertEqual(gf.qualifiers["ref"], "A")
        self.assertEqual(gf.qualifiers["source"], "ClinVar")
        self.assertEqual(gf.qualifiers["uid"], 246854)
        session.close()
        self.engine.dispose()


if __name__ == '__main__':
    unittest.main()
