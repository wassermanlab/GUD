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
from .genomicFeatureMixin1 import GFMixin1
from sqlalchemy.ext.declarative import declared_attr

class CNV(GFMixin1, Base):
    __tablename__ = "copy_number_variants"

    copy_number_change  = Column("copy_number_change", Integer, nullable=False)
    clinical_assertion  = Column("clinical_assertion", mysql.LONGBLOB, nullable=False)
    clinvar_accession   = Column("clinvar_accession", mysql.LONGBLOB, nullable=False)
    dbVar_accession     = Column("dbVar_accession", mysql.LONGBLOB, nullable=False)
    
    @declared_attr
    def __table_args__(cls):
        return (
        Index("ix_source_id", cls.source_id),
        Index("ix_cnv_region_id", cls.region_id),
        Index("ix_cnv_uid", cls.uid),

        {
            "mysql_engine": "MyISAM",
            "mysql_charset": "utf8"
        }
    )

    #not included in REST API
    @classmethod
    def is_unique(cls, session, regionID, sourceID, copy_number_change):
        q = session.query(cls).filter(cls.region_id == regionID,
                                      cls.source_id == sourceID, 
                                      cls.copy_number_change == copy_number_change)
        q = q.all()
        return len(q) == 0

    @classmethod
    def as_genomic_feature(self, feat):

        qualifiers = {   
            "uid": feat.CNV.uid, 
            "source": feat.Source.name, 
            "copy_number": feat.CNV.copy_number, 
            "clinical_interpretation": feat.CNV.clinical_interpretation, 
            "variant_type": feat.CNV.variant_type    
        }
        return GenomicFeature(
            feat.Region.chrom,
            int(feat.Region.start) ,
            int(feat.Region.end),
            strand = feat.Region.strand,
            feat_type = "CopyNumberVariant",
            feat_id = "%s_%s"%(self.__tablename__, feat.CNV.uid),
            qualifiers = qualifiers)
