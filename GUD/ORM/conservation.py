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

class Conservation(Base):

    __tablename__ = "conservation"

    uid = Column(
        "uid",
        mysql.INTEGER(unsigned=True)
    )

    regionID = Column(
        "regionID",
        Integer,
        ForeignKey("regions.uid"),
        nullable=False
    )

    score = Column("score", Float)

    sourceID = Column(
        "sourceID",
        Integer,
        ForeignKey("sources.uid"),
        nullable=False
    )

    __table_args__ = (
        PrimaryKeyConstraint(uid),
        UniqueConstraint(
            regionID,
            sourceID
        ),
        Index("ix_regionID", regionID), # query by bin range
        {
            "mysql_engine": "MyISAM",
            "mysql_charset": "utf8"
        }
    )

    @classmethod
    def is_unique(cls, session, regionID,
        sourceID):

        q = session.query(cls)\
            .filter(
                cls.regionID == regionID,
                cls.sourceID == sourceID
            )

        return len(q.all()) == 0

    @classmethod
    def select_unique(cls, session, regionID,
        sourceID):

        q = session.query(cls)\
            .filter(
                cls.regionID == regionID,
                cls.sourceID == sourceID
            )

        return q.first()

    @classmethod
    def select_by_location(cls, session, chrom,
        start, end, as_genomic_feature=False):
        """
        Query objects by genomic location.
        """

        bins = Region._compute_bins(start, end)

        q = session.query(cls, Region, Source)\
            .join()\
            .filter(
                Region.uid == cls.regionID,
                Source.uid == cls.sourceID,
            )\
            .filter(
                Region.chrom == chrom,
                Region.start < end,
                Region.end > start
            )\
            .filter(Region.bin.in_(bins))

        if as_genomic_feature:

            feats = []

            # For each feature...
            for feat in q.all():
                feats.append(
                    cls.__as_genomic_feature(feat)
                )

            return feats
    
        return q.all()

    @classmethod
    def __as_genomic_feature(self, feat):
        
        qualifiers = {
            "uid": feat.Conservation.uid,
            "regionID": feat.Conservation.regionID, 
            "score": feat.Conservation.score,
            "sourceID": feat.Conservation.sourceID, 
            }

        return GenomicFeature(
            feat.Region.chrom,
            int(feat.Region.start),
            int(feat.Region.end),
            feat_type = self.__tablename__,
            feat_id = "%s_%s" % \
                (
                    self.__tablename__,
                    feat.Conservation.uid
                ),
            qualifiers = qualifiers
        )

    def __repr__(self):

        return "<%s(%s, %s, %s, %s)>" % \
            (
                self.__tablename__,
                "uid={}".format(self.uid),
                "regionID={}".format(self.regionID),
                "score={}".format(self.score),
                "sourceID={}".format(self.sourceID)
            )