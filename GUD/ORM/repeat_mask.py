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

class RepeatMask(Base):

    __tablename__ = "rmsk"

    uid = Column("uid", mysql.INTEGER(unsigned=True))
    regionID = Column("regionID", Integer, ForeignKey('regions.uid'), nullable=False)
    sourceID = Column("sourceID", Integer, ForeignKey('sources.uid'), nullable=False)
    swScore = Column("swScore", Float, nullable=False)
    repName = Column("repName", String(75), nullable=False)
    repClass = Column("repClass", String(75), nullable=False)
    repFamily = Column("repFamily", String(75), nullable=False)
    strand = Column("strand", mysql.CHAR(1), nullable=False)

    __table_args__ = (
        PrimaryKeyConstraint(uid),
        UniqueConstraint(regionID, sourceID),

        Index("ix_rmsk", regionID),

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
        q = session.query(cls).filter(cls.regionID == regionID, cls.sourceID == sourceID)
        q = q.all()
        if len(q) == 0:
            return True
        else: 
            return False 

    def __str__(self):
        return "{}\t{}\t{}\t{}".format(self.swScore, self.repName, self.repClass, self.strand)

    def __repr__(self):
        return "<RMSK(uid={}, regionID={}, sourceID={}, swScore={}, repName={}, repClass={}, strand={})>".format(
            self.uid, self.regionID, self.sourceID, self.swScore, self.repName, self.repClass, self.strand)