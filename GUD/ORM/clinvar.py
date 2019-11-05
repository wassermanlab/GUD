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
from .genomicFeatureMixin1 import GFMixin1
from sqlalchemy.ext.declarative import declared_attr

class ClinVar(GFMixin1, Base):

    __tablename__ = "clinvar"

    ##fields
    ref = Column("ref", mysql.VARCHAR(5000), nullable=False) 
    alt = Column("alt", mysql.VARCHAR(5000), nullable=False) 
    clinvar_variation_ID = Column("clinvar_variation_ID", mysql.INTEGER(unsigned=True), nullable=False)
    ##info
    ANN_Annotation = Column("ANN_Annotation", mysql.VARCHAR(250))
    ANN_Annotation_Impact = Column("ANN_Annotation_Impact", mysql.VARCHAR(50))
    ANN_Gene_Name = Column("ANN_Gene_Name", mysql.VARCHAR(50))
    ANN_Gene_ID = Column("ANN_Gene_ID", mysql.VARCHAR(100))
    ANN_Feature_Type = Column("ANN_Feature_Type", mysql.VARCHAR(500))
    ANN_Feature_ID = Column("ANN_Feature_ID", mysql.VARCHAR(50))
    CADD = Column("CADD", Float)
    CLNDISDB = Column("CLNDISDB", mysql.VARCHAR(3000))  
    CLNDN = Column("CLNDN", mysql.VARCHAR(3000))
    CLNSIG = Column("CLNSIG", mysql.VARCHAR(500))
    gnomad_exome_af_global = Column("gnomad_exome_af_global", Float)
    gnomad_exome_hom_global = Column("gnomad_exome_hom_global", Float)
    gnomad_genome_af_global = Column("gnomad_genome_af_global", Float)
    gnomad_genome_hom_global = Column("gnomad_genome_hom_global", Float) 

    @declared_attr
    def __table_args__(cls):
        return (
        UniqueConstraint(cls.clinvar_variation_ID),
        Index("ix_join", cls.source_id, cls.region_id),
        Index("ix_clinvar_id", cls.clinvar_variation_ID),
        {"mysql_engine": "InnoDB", "mysql_charset": "utf8"}
    )
   
    @classmethod
    def select_by_clinvarID(cls, session, query, clinvarIDs):
        """
        filter query by location
        """
        if (query is None):
            q = cls.make_query(session)
        else: 
            q = query 
        q = q.filter(cls.clinvar_variation_ID.in_(clinvarIDs))
        return q

    # not included in REST
    @classmethod
    def is_unique(cls, session, clinvarID):
        q = session.query(cls).filter(cls.clinvar_variation_ID == clinvarID)
        return len(q.all()) == 0

    @classmethod
    def as_genomic_feature(self, feat):
        # Define qualifiers
        qualifiers = {
            "uid": feat.ClinVar.uid,
            "source": feat.Source.name,
            "ref": feat.ClinVar.ref,
            "alt": feat.ClinVar.alt,
            "clinvarID": feat.ClinVar.clinvar_variation_ID,
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
            feat_id = "%s_%s"%(self.__tablename__, feat.ClinVar.uid),
            qualifiers = qualifiers)

 