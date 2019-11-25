from binning import containing_bins, contained_bins
from sqlalchemy import Column, Index, PrimaryKeyConstraint, String, ForeignKey, UniqueConstraint, Integer, Float
from sqlalchemy.dialects import mysql
from .genomicFeatureMixin1 import GFMixin1
from .base import Base
from .region import Region
from .source import Source
from sqlalchemy.ext.declarative import declared_attr

class Mask(GFMixin1, Base):
    # table declerations
    __tablename__ = "masks"
    
    name = Column("name", String(75))
    score = Column("score", Float)

    @declared_attr
    def __table_args__(cls):
        return (
            UniqueConstraint(cls.region_id, cls.source_id),
            Index("ix_join", cls.source_id, cls.region_id),
            Index("ix_clinvar_id", cls.name),
            {"mysql_engine": "InnoDB", "mysql_charset": "utf8"}
        )

    # class methods
    @classmethod
    def as_genomic_feature(self, feat):
        """
        extend parent class by adding qualifiers
        """
        qualifiers = {
            "uid": feat.Mask.uid,
            "name": feat.Mask.name,
            "source": feat.Source.name,
        }
        genomic_feature = super().as_genomic_feature(feat)
        genomic_feature.qualifiers = qualifiers
        genomic_feature.score = feat.Mask.score
        return genomic_feature