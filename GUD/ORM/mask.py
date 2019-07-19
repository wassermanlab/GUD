from binning import containing_bins, contained_bins
from sqlalchemy import Column, Index, PrimaryKeyConstraint, String, ForeignKey, UniqueConstraint, Integer, Float
from sqlalchemy.dialects import mysql

from .base import Base
from .region import Region
from .source import Source

class Mask(Base):

    __tablename__ = "masks"

    uid = Column("uid", mysql.INTEGER(unsigned=True))
    regionID = Column("regionID", Integer, ForeignKey("regions.uid"), nullable=False)
    sourceID = Column("sourceID", Integer, ForeignKey("sources.uid"), nullable=False)
    name = Column("name", String(75))
    score = Column("score", Float)

    __table_args__ = (
        PrimaryKeyConstraint(uid),
        UniqueConstraint(regionID, sourceID),
        Index("ix_regionID", regionID), # query by bin range
        {
            "mysql_engine": "MyISAM",
            "mysql_charset": "utf8"
        }
    )

    @classmethod
    def is_unique(cls, session, regionID, sourceID):

        q = session.query(cls).filter(cls.regionID == regionID, cls.sourceID == sourceID)

        return len(q.all()) == 0

    @classmethod
    def select_unique(cls, session, regionID, sourceID):

        q = session.query(cls).filter(cls.regionID == regionID, cls.sourceID == sourceID)

        return q.first()

    def __repr__(self):

        return "<%s(%s, %s, %s, %s, %s)>" % \
            (   
                self.__tablename__,
                "uid={}".format(self.uid),
                "regionID={}".format(self.regionID),
                "sourceID={}".format(self.sourceID),
                "name={}".format(self.name),
                "score={}".format(self.score),
            )