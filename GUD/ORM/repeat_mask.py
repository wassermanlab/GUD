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


class RepeatMask(GFMixin1, Base):
    # table declerations 
    __tablename__ = "rmsk"

    swScore = Column("swScore", Float, nullable=False)
    repName = Column("repName", String(75), nullable=False)
    repClass = Column("repClass", String(75), nullable=False)
    repFamily = Column("repFamily", String(75), nullable=False)
    
    @declared_attr
    def __table_args__(cls):
        return (
        UniqueConstraint(cls.region_id, cls.source_id),
        Index("ix_join", cls.region_id, cls.source_id),
        {"mysql_engine": "InnoDB", "mysql_charset": "utf8"}
    )

    # class methods 
    @classmethod
    def as_genomic_feature(self, feat):
        """
        extend parent class by adding qualifiers
        """
        qualifiers = {
            "uid": feat.RepeatMask.uid,
            "source": feat.Source.name,
            "swScore": feat.RepeatMask.swScore,
            "repName": feat.RepeatMask.repName,
            "repClass": feat.RepeatMask.repClass,
            "repFamily": feat.RepeatMask.repFamily,
            "repClass": feat.RepeatMask.repClass
        }
        genomic_feature = super().as_genomic_feature(feat)
        genomic_feature.qualifiers = qualifiers
        return genomic_feature
