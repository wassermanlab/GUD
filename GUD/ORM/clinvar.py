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

    @declared_attr
    def __table_args__(cls):
        return (
        UniqueConstraint(cls.clinvarID),

        Index("ix_clinvar", cls.region_id),
        Index("ix_clinvar_id", cls.clinvarID),

        {  
            "mysql_engine": "MyISAM",
            "mysql_charset": "utf8"
        }
    )
   
    @classmethod
    def select_by_clinvarID(cls, session, clinvarID, limit, offset):
        """
        Query clinvar objects by clinvarID. If no name is provided,
        query all ids.
        """
        q = session.query(cls)
        q = session.query(cls, Region, Source).\
        join().\
        filter(Region.uid == cls.region_id,
                Source.uid == cls.source_id,
                cls.clinvarID == clinvarID)
        
        return (q.count(), q.offset(offset).limit(limit))

    # not included in REST
    @classmethod
    def is_unique(cls, session, clinvarID):
        q = session.query(cls).filter(cls.clinvarID == clinvarID)
        return len(q.all()) == 0

    @classmethod
    def as_genomic_feature(self, feat):
        # Define qualifiers
        qualifiers = {
            "uid": feat.ClinVar.uid,
            "source": feat.Source.name,
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
            feat_id = "%s_%s"%(self.__tablename__, feat.ClinVar.uid),
            qualifiers = qualifiers)

 