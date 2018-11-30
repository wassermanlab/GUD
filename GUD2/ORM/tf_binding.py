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

class TFBinding(Base):

    __tablename__ = "tf_binding"

    uid = Column("uid", mysql.INTEGER(unsigned=True))
    regionID = Column("regionID", Integer, ForeignKey('regions.uid'), nullable=False)
    sourceID = Column("sourceID", Integer, ForeignKey('sources.uid'), nullable=False)
    sampleID = Column("sampleID", Integer, ForeignKey('samples.uid'), nullable=False)
    experimentID = Column("experimentID", Integer, ForeignKey('experiments.uid'), nullable=False)
    tf = Column("tf", String(250), nullable=False)

    __table_args__ = (
        PrimaryKeyConstraint(uid),
        UniqueConstraint(regionID, sourceID, sampleID, experimentID, tf),

        Index("ix_tf_binding", regionID), 

        {
            "mysql_engine": "MyISAM",
            "mysql_charset": "utf8"
        }
    )

    @classmethod
    def select_unique(cls, session, regionID, sourceID, sampleID, experimentID, tf):
        """
        Query objects by name of sample type. 
        """
        q = session.query(cls).filter(
            cls.regionID == regionID, 
            cls.sourceID == sourceID,
            cls.sampleID == sampleID,
            cls.experimentID == experimentID, 
            cls.tf == tf)

        return q.first()

    def __str__(self):
        return "{}".format(self.tf)

    def __repr__(self):
        return "<TFBinding(uid={}, regionID={}, sourceID={}, sampleID={}, experimentID={}, tf={})>".format(
            self.uid, self.regionID, self.sourceID, self.sampleID, self.experimentID, self.tf)