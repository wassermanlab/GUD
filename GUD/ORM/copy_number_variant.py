from binning import (
    containing_bins,
    contained_bins,
    assign_bin
)
from sqlalchemy import (
    Column, Index, PrimaryKeyConstraint, String, ForeignKey,
    UniqueConstraint, CheckConstraint, Integer
)
from sqlalchemy.dialects import mysql

from .base import Base
from .region import Region
from .source import Source
from .genomic_feature import GenomicFeature
from .genomicFeature import GF
from sqlalchemy.ext.declarative import declared_attr

class CNV(GF, Base):
    __tablename__ = "copy_number_variants"

    variant_type = Column("variant_type", String(50), nullable=False)
    copy_number = Column("copy_number", Integer, nullable=False)
    clinical_interpretation = Column("clinical_interpretation", String(50), nullable=False)


    @declared_attr
    def __table_args__(cls):
        return (
        Index("ix_cnv", cls.region_id),
        Index("ix_cnv_uid", cls.uid),
        Index("ix_cnv_clinical_interpretation", cls.clinical_interpretation),

        {
            "mysql_engine": "MyISAM",
            "mysql_charset": "utf8"
        }
    )

    #not included in REST API
    @classmethod
    def is_unique(cls, session, name):
        q = session.query(cls).filter(cls.uid == name)
        q = q.all()
        return len(q) == 0

    @classmethod
    def as_genomic_feature(self, feat):

        qualifiers = {   
            "uid": feat.CNV.uid, 
            "regionID": feat.CNV.region_id, 
            "sourceID": feat.CNV.source_id, 
            "copy_number": feat.CNV.copy_number, 
            "clinical_interpretation": feat.CNV.clinical_interpretation, 
            "variant_type": feat.CNV.variant_type    
        }
        return GenomicFeature(
            feat.Region.chrom,
            int(feat.Region.start) + 1,
            int(feat.Region.end),
            strand = feat.Region.strand,
            feat_type = "CopyNumberVariant",
            feat_id = "%s_%s"%(self.__tablename__, feat.CNV.uid),
            qualifiers = qualifiers)
