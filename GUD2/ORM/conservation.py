from sqlalchemy import (
    Column, Index, PrimaryKeyConstraint, String, ForeignKey,
    UniqueConstraint, Integer, Float
)
from sqlalchemy.dialects import mysql
from GUD2.ORM.region import Region
from GUD2.ORM.source import Source
from GUD2.ORM.base import Base
from binning import containing_bins, contained_bins

class Conservation(Base):

    __tablename__ = "conservation"

    uid = Column("uid", mysql.INTEGER(unsigned=True))
    regionID = Column("regionID", Integer, ForeignKey("regions.uid"), nullable=False)
    sourceID = Column("sourceID", Integer, ForeignKey("sources.uid"), nullable=False)
    score = Column("score", Float)

    __table_args__ = (
        PrimaryKeyConstraint(uid),
        UniqueConstraint(regionID, sourceID),

        Index("ix_cons", regionID),
        Index("ix_cons_score", score),

        {
            "mysql_engine": "MyISAM",
            "mysql_charset": "utf8"
        }
    )

    @classmethod
    def select_by_location(cls, session, chrom, start, end):
        """
        Query objects based off of their location being within the start only
        motifs through that  
         """

        bins = set(containing_bins(start, end) + contained_bins(start, end))
        
        q = session.query(cls, Region).\
        join().\
        filter(Region.uid == cls.regionID).\
        filter(Region.chrom == chrom, Region.end > start, Region.start < end).\
        filter(Region.bin.in_(bins))
        return q.all()

    @classmethod
    def is_unique(cls, session, regionID, sourceID):

        q = session.query(cls).filter(
            cls.regionID == regionID,
            cls.sourceID == sourceID
        )

        return len(q.all()) == 0

    def __str__(self):
        return "{}".format(self.score)

    def __repr__(self):
        return "<Conservation(uid={}, regionID={}, sourceID={}, score={})>".format(
            self.uid, self.regionID, self.sourceID, self.score)