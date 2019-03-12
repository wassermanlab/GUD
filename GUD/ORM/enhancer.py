from binning import (
    containing_bins,
    contained_bins
)
from sqlalchemy import (
    Column, Index, PrimaryKeyConstraint, ForeignKey,
    UniqueConstraint, Integer
)
from sqlalchemy.dialects import mysql

from .base import Base
from .experiment import Experiment
from .region import Region
from .sample import Sample
from .source import Source

class Enhancer(Base):

    __tablename__ = "enhancers"

    uid = Column("uid", mysql.INTEGER(unsigned=True))
    regionID = Column("regionID", Integer, ForeignKey("regions.uid"), nullable=False)
    sourceID = Column("sourceID", Integer, ForeignKey("sources.uid"), nullable=False)
    sampleID = Column("sampleID", Integer, ForeignKey("samples.uid"), nullable=False)
    experimentID = Column("experimentID", Integer, ForeignKey("experiments.uid"), nullable=False) 

    __table_args__ = (
        PrimaryKeyConstraint(uid),
        UniqueConstraint(regionID, sourceID, sampleID, experimentID),

        Index("ix_enhancer", regionID),
        Index("ix_enhancer_sample", sampleID),

        {
            "mysql_engine": "MyISAM",
            "mysql_charset": "utf8"
        }
    )

    @classmethod
    def is_unique(cls, session, regionID, sourceID, sampleID, experimentID):

        q = session.query(cls).filter(
            cls.regionID == regionID,
            cls.sourceID == sourceID,
            cls.sampleID == sampleID,
            cls.experimentID == experimentID
        )

        return len(q.all()) == 0

#    @classmethod
#    def select_by_location(cls, session, chrom, start, end, samples=[]):
#        """
#        Query objects based off of their location being within
#        the start only motifs through that  
#        """
#
#        bins = set(containing_bins(start, end) + contained_bins(start, end))
#
#        q = session.query(cls, Region, Sample).\
#        join().\
#        filter(Region.uid == cls.regionID).\
#        filter(Region.chrom == chrom, Region.end > start, Region.start < end).\
#        filter(Region.bin.in_(bins))
#
#        if samples:
#            q = q.filter()
#
#        return q.all()

#    @classmethod
#    def select_unique(cls, session, regionID, sourceID, sampleID, experimentID):
#        """
#        Query objects by name of sample type. 
#        """
#
#        q = session.query(cls).filter(
#            cls.regionID == regionID,
#            cls.sourceID == sourceID,
#            cls.sampleID == sampleID,
#            cls.experimentID == experimentID
#        )
#
#        return q.first()

    def __repr__(self):
        return "<Enhancer(uid={}, regionID={}, sourceID={}, sampleID={}, experimentID={})>".format(
            self.uid, self.regionID, self.sourceID, self.sampleID, self.experimentID)