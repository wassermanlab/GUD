from sqlalchemy import (
    Column, Index, PrimaryKeyConstraint, String, ForeignKey,
    UniqueConstraint, CheckConstraint, Integer, Float
)
from sqlalchemy.dialects import mysql
from GUD2.ORM.region import Region
from GUD2.ORM.source import Source
from GUD2.ORM.sample import Sample
from GUD2.ORM.experiment import Experiment
from GUD2.ORM.base import Base
from binning import containing_bins, contained_bins

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

        Index("ix_enhancer", regionID), ## query by bin range 
        Index("ix_enhancer_sample", sampleID),

        {
            "mysql_engine": "MyISAM",
            "mysql_charset": "utf8"
        }
    )

    @classmethod
    def select_unique(cls, session, regionID, sourceID, sampleID, experimentID):
        """
        Query objects by name of sample type. 
        """

        q = session.query(cls).filter(
            cls.regionID == regionID,
            cls.sourceID == sourceID,
            cls.sampleID == sampleID,
            cls.experimentID == experimentID
        )

        return q.first()

    @classmethod
    def is_unique(cls, session, regionID, sourceID, sampleID, experimentID):

        q = session.query(cls).filter(
            cls.regionID == regionID,
            cls.sourceID == sourceID,
            cls.sampleID == sampleID,
            cls.experimentID == experimentID
        )

        return len(q.all()) == 0

    def __repr__(self):
        return "<Enhancer(uid={}, regionID={}, sourceID={}, sampleID={}, experimentID={})>".format(
            self.uid, self.regionID, self.sourceID, self.sampleID, self.experimentID)