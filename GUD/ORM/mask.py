from binning import containing_bins, contained_bins
from sqlalchemy import Column, Index, PrimaryKeyConstraint, String, ForeignKey, UniqueConstraint, Integer, Float
from sqlalchemy.dialects import mysql
from .genomicFeatureMixin1 import GFMixin1
from .base import Base
from .region import Region
from .source import Source

class Mask(GFMixin1 Base):
    # table declerations
    __tablename__ = "masks"
    
    name = Column("name", String(75))
    score = Column("score", Float)

    __table_args__ = (
        Index("ix_regionID", regionID), # query by bin range
        {"mysql_engine": "InnoDB","mysql_charset": "utf8"}
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