from binning import (
    containing_bins,
    contained_bins
)
from sqlalchemy import (
    Column, Index, PrimaryKeyConstraint, String, ForeignKey,
    UniqueConstraint, CheckConstraint, Integer, Float
)
from sqlalchemy.dialects import mysql

from .base import Base
from .region import Region
from .source import Source
from .genomic_feature import GenomicFeature

class ClinVar(Base):

    __tablename__ = "clinvar"

    uid = Column("uid", mysql.INTEGER(unsigned=True))
    regionID = Column("regionID", Integer, ForeignKey('regions.uid'), nullable=False)
    sourceID = Column("sourceID", Integer, ForeignKey('sources.uid'), nullable=False)
    ##fields
    ref = Column("ref", mysql.VARCHAR(5000), nullable=False) ## check these
    alt = Column("alt", mysql.VARCHAR(5000), nullable=False) ## check these
    clinvarID = Column("clinvarID", mysql.VARCHAR(7), nullable=False)
    ##info
    ANN_Annotation = Column("ANN_Annotation", mysql.VARCHAR(250))
    ANN_Annotation_Impact = Column("ANN_Annotation_Impact", mysql.VARCHAR(50))
    ANN_Gene_Name = Column("ANN_Gene_Name", mysql.VARCHAR(50))
    ANN_Gene_ID = Column("ANN_Gene_ID", mysql.VARCHAR(50))
    ANN_Feature_Type = Column("ANN_Feature_Type", mysql.VARCHAR(500))
    ANN_Feature_ID = Column("ANN_Feature_ID", mysql.VARCHAR(50))
    CADD = Column("CADD", Float)
    CLNDISDB = Column("CLNDISDB", mysql.VARCHAR(3000)) ## check 
    CLNDN = Column("CLNDN", mysql.VARCHAR(3000))
    CLNSIG = Column("CLNSIG", mysql.VARCHAR(500))
    gnomad_exome_af_global = Column("gnomad_exome_af_global", Float)
    gnomad_exome_hom_global = Column("gnomad_exome_hom_global", Float)
    gnomad_genome_af_global = Column("gnomad_genome_af_global", Float)
    gnomad_genome_hom_global = Column("gnomad_genome_hom_global", Float) 

    __table_args__ = (
        PrimaryKeyConstraint(uid),
        UniqueConstraint(clinvarID),

        Index("ix_clinvar", regionID),
        Index("ix_clinvar_id", clinvarID),

        {  
            "mysql_engine": "MyISAM",
            "mysql_charset": "utf8"
        }
    )

    @classmethod
    def select_by_location(cls, session, chrom,
        start, end, as_genomic_feature=False):
        """
        Query objects by genomic location.
        """

        bins = Region._compute_bins(start, end)

        q = session.query(cls, Region, Source)\
            .join()\
            .filter(
                Region.uid == cls.regionID,
                Source.uid == cls.sourceID,
            )\
            .filter(
                Region.chrom == chrom,
                Region.start < end,
                Region.end > start
            )\
            .filter(Region.bin.in_(bins))

        if as_genomic_feature:

            feats = []

            # For each feature...
            for feat in q.all():
                feats.append(
                    cls.__as_genomic_feature(feat)
                )
            return feats
        return q.all()
   
    @classmethod
    def select_by_name(cls, session, clinvarID, as_genomic_feature=False):
        """
        Query refGene objects by common name. If no name is provided,
        query all genes.
        """
        q = session.query(cls)
        q = session.query(cls, Region).\
        join().\
        filter(Region.uid == cls.regionID).\
        filter(cls.clinvarID == clinvarID)
        
        if as_genomic_feature:
            return cls.__as_genomic_feature(q.first())
        return q.first()

    @classmethod
    def is_unique(cls, session, clinvarID):
        q = session.query(cls).filter(cls.clinvarID == clinvarID)
        return len(q.all()) == 0

    @classmethod
    def __as_genomic_feature(self, feat):
        # Define qualifiers
        qualifiers = {
            "ref": feat.ClinVar.ref,
            "alt": feat.ClinVar.alt,
            "clinvarID": feat.ClinVar.clinvarID,
            "ANN_Annotation": feat.ClinVar.ANN_Annotation,
            "ANN_Annotation_Impact": feat.ClinVar.ANN_Annotation_Impact,
            "ANN_Gene_Name": feat.ClinVar.ANN_Gene_Name,
            "ANN_Gene_ID": feat.ClinVar.ANN_Gene_ID,
            "ANN_Feature_Type": feat.ClinVar.ANN_Feature_Type,
            "ANN_Feature_ID": feat.ClinVar.ANN_Feature_ID,
            "CADD": feat.ClinVar.CADD,
            "CLNDISDB": feat.ClinVar.CLNDISDB,
            "CLNDN": feat.ClinVar.CLNDN,
            "CLNSIG": feat.ClinVar.CLNSIG,
            "gnomad_exome_af_global": feat.ClinVar.gnomad_exome_af_global,
            "gnomad_exome_hom_global": feat.ClinVar.gnomad_exome_hom_global,
            "gnomad_genome_af_global": feat.ClinVar.gnomad_genome_af_global,
            "gnomad_genome_hom_global": feat.ClinVar.gnomad_genome_hom_global,          
        }
        return GenomicFeature(
            feat.Region.chrom,
            int(feat.Region.start),
            int(feat.Region.end),
            strand = feat.Region.strand,
            feat_type = "ClinVar",
            feat_id = feat.ClinVar.clinvarID, 
            qualifiers = qualifiers)

    def __repr__(self):
        return "<ClinVar(uid={}, regionID={}, sourceID={},\
        ref={}, alt={}, clinvarID={},\
        ANN_Annotation={}, ANN_Annotation_Impact={}, ANN_Gene_Name={},\
        ANN_Gene_ID={}, ANN_Feature_Type={}, ANN_Feature_ID={},\
        CADD={}, CLNDISDB={}, CLNDN={}, CLNSIG={},\
        gnomad_exome_af_global={}, gnomad_exome_hom_global={}, gnomad_genome_af_global={}, gnomad_genome_hom_global={})>".format(
            self.uid, self.regionID, self.sourceID, self.ref, self.alt, self.clinvarID,
            self.ANN_Annotation, self.ANN_Annotation_Impact, self.ANN_Gene_Name,
            self.ANN_Gene_ID, self.ANN_Feature_Type, self.ANN_Feature_ID,
            self.CADD, self.CLNDISDB, self.CLNDN, self.CLNSIG, 
            self.gnomad_exome_af_global, self.gnomad_exome_hom_global,
            self.gnomad_genome_af_global, self.gnomad_genome_hom_global)
