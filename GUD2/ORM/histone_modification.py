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

class HistoneModification(Base):

    __tablename__ = "histone_modification"

    uid = Column("uid", mysql.INTEGER(unsigned=True))
    regionID = Column("regionID", Integer, ForeignKey('regions.uid'), nullable=False)
    sourceID = Column("sourceID", Integer, ForeignKey('sources.uid'), nullable=False)
    sampleID = Column("sampleID", Integer, ForeignKey('samples.uid'), nullable=False)
    experimentID = Column("experimentID", Integer, ForeignKey('experiments.uid'), nullable=False)
    histone_type = Column("histone_type", String(25), nullable=False)

    __table_args__ = (
        PrimaryKeyConstraint(uid),
        UniqueConstraint(regionID, sourceID, sampleID, experimentID, histone_type),

        Index("ix_tf_binding", regionID),

        {
            "mysql_engine": "MyISAM",
            "mysql_charset": "utf8"
        }
    )

    @classmethod
    def select_unique(cls, session, regionID, sourceID, sampleID, experimentID, histone_type):
        """
        Query objects by name of sample type. 
        """
        q = session.query(cls).filter(
            cls.regionID == regionID, 
            cls.sourceID == sourceID,
            cls.sampleID == sampleID,
            cls.experimentID == experimentID, 
            cls.histone_type == histone_type)

        return q.first()

    def __str__(self):
        return "{}".format(self.histone_type)

    def __repr__(self):
        return "<HistoneModification(uid={}, regionID={}, sourceID={}, sampleID={}, experimentID={}, histone_type={})>".format(
            self.uid, self.regionID, self.sourceID, self.sampleID, self.experimentID, self.histone_type)