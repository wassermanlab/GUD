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
from .sample import Sample
from .experiment import Experiment

class DNAAccessibility(Base):

    __tablename__ = "dna_accessibility"

    uid = Column("uid", mysql.INTEGER(unsigned=True))
    regionID = Column("regionID", Integer, ForeignKey("regions.uid"), nullable=False)
    sourceID = Column("sourceID", Integer, ForeignKey("sources.uid"), nullable=False)
    sampleID = Column("sampleID", Integer, ForeignKey("samples.uid"), nullable=False)
    experimentID = Column("experimentID", Integer, ForeignKey("experiments.uid"), nullable=False)

    __table_args__ = (
        PrimaryKeyConstraint(uid),
        UniqueConstraint(regionID, sourceID, sampleID, experimentID),

        Index("ix_dna_accessibility", regionID),

        {
            "mysql_engine": "MyISAM",
            "mysql_charset": "utf8"
        }
    )

    @classmethod
    def is_unique(cls, session, regionID, sourceID, sampleID,
        experimentID):

        q = session.query(cls).filter(
            cls.regionID == regionID,
            cls.sourceID == sourceID,
            cls.sampleID == sampleID,
            cls.experimentID == experimentID
        )

        return len(q.all()) == 0

    @classmethod
    def select_unique(cls, session, regionID, sourceID, sampleID,
        experimentID):

        q = session.query(cls).filter(
            cls.regionID == regionID, 
            cls.sourceID == sourceID,
            cls.sampleID == sampleID,
            cls.experimentID == experimentID)

        return q.first()

    def __repr__(self):
        return "<DNAAccessibility(uid={}, regionID={}, sourceID={}, sampleID={}, experimentID={})>".format(
            self.uid, self.regionID, self.sourceID, self.sampleID, self.experimentID)