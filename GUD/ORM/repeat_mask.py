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

    __tablename__ = "rmsk"

    swScore = Column("swScore", Float, nullable=False)
    repName = Column("repName", String(75), nullable=False)
    repClass = Column("repClass", String(75), nullable=False)
    repFamily = Column("repFamily", String(75), nullable=False)
    __table_args__ = (
        UniqueConstraint(cls.region_id, cls.source_id),
        Index("ix_join", region_id, source_id),
        {"mysql_engine": "InnoDB", "mysql_charset": "utf8"}
    )

    @classmethod
    def is_unique(cls, session, regionID, sourceID):
        q = session.query(cls).filter(
            cls.region_id == regionID, cls.source_id == sourceID)
        q = q.all()
        if len(q) == 0:
            return True
        else:
            return False

    @classmethod
    def as_genomic_feature(self, feat):
        # Define qualifiers
        qualifiers = {
            "uid": feat.RepeatMask.uid,
            "source": feat.Source.name,
            "swScore": feat.RepeatMask.swScore,
            "repName": feat.RepeatMask.repName,
            "repClass": feat.RepeatMask.repClass,
            "repFamily": feat.RepeatMask.repFamily,
            "repClass": feat.RepeatMask.repClass
        }
        return GenomicFeature(
            feat.Region.chrom,
            int(feat.Region.start),
            int(feat.Region.end),
            strand=feat.Region.strand,
            feat_type="RepeatMask",
            feat_id="%s_%s" % (self.__tablename__, feat.RepeatMask.uid),
            qualifiers=qualifiers)
