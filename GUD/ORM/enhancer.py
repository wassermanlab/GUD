from binning import (
    containing_bins,
    contained_bins
)
from sqlalchemy import (
    Column,
    Index,
    Integer,
    PrimaryKeyConstraint,
    ForeignKey,
    UniqueConstraint
)
from sqlalchemy.dialects import mysql

from .base import Base
from .experiment import Experiment
from .region import Region
from .sample import Sample
from .source import Source

class Enhancer(Base):

    __tablename__ = "enhancers"

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

    sampleID = Column(
        "sampleID",
        Integer,
        ForeignKey("samples.uid"),
        nullable=False
    )

    experimentID = Column(
        "experimentID",
        Integer,
        ForeignKey("experiments.uid"),
        nullable=False
    )

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
            sampleID,
            experimentID,
            sourceID
        ),
        Index("ix_regionID", regionID), # query by bin range
        Index("ix_sampleID", sampleID),
        {
            "mysql_engine": "MyISAM",
            "mysql_charset": "utf8"
        }
    )

    @classmethod
    def is_unique(cls, session, regionID, sampleID,
        experimentID, sourceID):

        q = session.query(cls).\
            filter(
                cls.regionID == regionID,
                cls.sampleID == sampleID,
                cls.experimentID == experimentID,
                cls.sourceID == sourceID
            )

        return len(q.all()) == 0

    @classmethod
    def select_unique(cls, session, regionID,
        sampleID, experimentID, sourceID):

        q = session.query(cls).\
            filter(
                cls.regionID == regionID,
                cls.sampleID == sampleID,
                cls.experimentID == experimentID,
                cls.sourceID == sourceID
            )

        return q.first()

    def __repr__(self):

        return "<Enhancer(%s, %s, %s, %s, %s)>" % \
            (
                "uid={}".format(self.uid),
                "regionID={}".format(self.regionID),
                "sampleID={}".format(self.sampleID),
                "experimentID={}".format(self.experimentID),
                "sourceID={}".format(self.sourceID)
            )