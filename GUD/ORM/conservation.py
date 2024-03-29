from sqlalchemy import (Column, Index, UniqueConstraint)
from sqlalchemy.dialects import mysql
from .base import Base
from .source import Source
from .genomicFeatureMixin1 import GFMixin1
from sqlalchemy.ext.declarative import declared_attr


class Conservation(GFMixin1, Base):
    # table declerations
    __tablename__ = "conservation"
    score = Column("score", mysql.FLOAT)

    @declared_attr
    def __table_args__(cls):
        return (
            UniqueConstraint(cls.region_id, cls.source_id),
            # query by bin range
            Index("ix_join", cls.source_id, cls.region_id),
            {"mysql_engine": "InnoDB", "mysql_charset": "utf8"}
        )
    # class methods
    @classmethod
    def as_genomic_feature(self, feat):
        """
        extend parent class by adding qualifiers
        """
        qualifiers = {
            "uid": feat.Conservation.uid,
            "source": feat.sourceName,
        }
        genomic_feature = super().as_genomic_feature(feat)
        genomic_feature.score = feat.Conservation.score
        genomic_feature.qualifiers = qualifiers
        return genomic_feature
