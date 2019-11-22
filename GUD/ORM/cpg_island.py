from sqlalchemy import (
    Column,
    Float,
    Index,
    Integer,
    UniqueConstraint
)

from .base import Base
from .region import Region
from .source import Source
from .genomic_feature import GenomicFeature
from .genomicFeatureMixin1 import GFMixin1
from sqlalchemy.ext.declarative import declared_attr

class CpGIsland(GFMixin1, Base):

    __tablename__ = "cpg_islands"

    cpgNum = Column("cpgNum", Integer)
    gcNum = Column("gcNum", Integer)
    perCpg = Column("perCpg", Float)
    perGc = Column("perGc", Float)
    obsExp = Column("obsExp", Float)

    @declared_attr
    def __table_args__(cls):
        return (
        UniqueConstraint(cls.region_id, cls.source_id),
        Index("ix_join", cls.source_id, cls.region_id), # query by bin range
        {"mysql_engine": "InnoDB", "mysql_charset": "utf8"}
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
            "uid": feat.CpGIsland.uid,
            "cpgNum": feat.CpGIsland.cpgNum,
            "gcNum": feat.CpGIsland.gcNum,
            "perCpg": feat.CpGIsland.perCpg,
            "perGc": feat.CpGIsland.perGc,
            "obsExp": feat.CpGIsland.obsExp,
            "source": feat.Source.name,
        }

        return GenomicFeature(
            feat.Region.chrom,
            int(feat.Region.start),
            int(feat.Region.end),
            feat_type="CpGIsland",
            feat_id="%s_%s" % (self.__tablename__, feat.CpGIsland.uid),
            qualifiers=qualifiers
        )