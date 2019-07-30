from sqlalchemy import (
    Column,
    Index,
    Float,
    ForeignKey,
    Integer,
    PrimaryKeyConstraint,
    UniqueConstraint
)
from sqlalchemy.dialects import mysql

from .base import Base
from .region import Region
from .source import Source
from .genomic_feature import GenomicFeature
from .genomicFeatureMixin1 import GFMixin1
from sqlalchemy.ext.declarative import declared_attr


class Conservation(GFMixin1, Base):

    __tablename__ = "conservation"

    score = Column("score", Float)

    @declared_attr
    def __table_args__(cls):
        return (
            UniqueConstraint(cls.region_id, cls.source_id),
            Index("ix_region_id", cls.region_id),  # query by bin range
            {
                "mysql_engine": "MyISAM",
                "mysql_charset": "utf8"
            }
        )

    @classmethod
    def is_unique(cls, session, regionID, sourceID):

        q = session.query(cls)\
            .filter(
                cls.region_id == regionID,
                cls.source_id == sourceID
        )

        return len(q.all()) == 0

    @classmethod
    def select_unique(cls, session, regionID, sourceID):

        q = session.query(cls)\
            .filter(
                cls.region_id == regionID,
                cls.source_id == sourceID
        )

        return q.first()

    @classmethod
    def as_genomic_feature(self, feat):

        qualifiers = {
            "uid": feat.Conservation.uid,
            "score": feat.Conservation.score,
            "source": feat.Source.name,
        }

        return GenomicFeature(
            feat.Region.chrom,
            int(feat.Region.start),
            int(feat.Region.end),
            feat_type="Conservation",
            feat_id="%s_%s" % (self.__tablename__, feat.Conservation.uid),
            qualifiers=qualifiers
        )