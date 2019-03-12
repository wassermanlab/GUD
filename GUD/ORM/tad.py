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
from .experiment import Experiment
from .region import Region
from .sample import Sample
from .source import Source

class TAD(Base):

    __tablename__ = "tad"

    uid = Column("uid", mysql.INTEGER(unsigned=True))
    regionID = Column("regionID", Integer, ForeignKey('regions.uid'), nullable=False)
    sourceID = Column("sourceID", Integer, ForeignKey('sources.uid'), nullable=False)
    sampleID = Column("sampleID", Integer, ForeignKey('samples.uid'), nullable=False)
    experimentID = Column("experimentID", Integer, ForeignKey('experiments.uid'), nullable=False)
    restriction_enzyme = Column("restriction_enzyme", String(25), nullable=False)

    __table_args__ = (
        PrimaryKeyConstraint(uid),
        UniqueConstraint(regionID, sourceID, sampleID, experimentID, restriction_enzyme),

        Index("ix_tad", regionID),

        {
            "mysql_engine": "MyISAM",
            "mysql_charset": "utf8"
        }
    )

    @classmethod
    def select_unique(cls, session, regionID, sourceID, sampleID, experimentID, restriction_enzyme):
        """
        Query objects by name of sample type. 
        """
        q = session.query(cls).filter(
            cls.regionID == regionID, 
            cls.sourceID == sourceID,
            cls.sampleID == sampleID,
            cls.experimentID == experimentID, 
            cls.restriction_enzyme == restriction_enzyme)

        return q.first()

    def __repr__(self):
        return "<TAD(uid={}, regionID={}, sourceID={}, sampleID={}, experimentID={}, restriction_enzyme={})>".format(
            self.uid, self.regionID, self.sourceID, self.sampleID, self.experimentID, self.restriction_enzyme)