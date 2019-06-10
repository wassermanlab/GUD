from binning import (
    containing_bins,
    contained_bins
)
from sqlalchemy import (
    Column,
    ForeignKey,
    Index,
    Integer,
    PrimaryKeyConstraint,
    String,
    UniqueConstraint
)
from sqlalchemy.dialects import mysql

from .base import Base
from .experiment import Experiment
from .region import Region
from .sample import Sample
from .source import Source

class HistoneModification(Base):

    __tablename__ = "histone_modifications"

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

    histone_type = Column(
        "histone_type",
        String(25),
        nullable=False
    )

    __table_args__ = (
        PrimaryKeyConstraint(uid),
        UniqueConstraint(
            regionID,
            sampleID,
            experimentID,
            sourceID,
            histone_type

        ),
        Index("ix_regionID", regionID), # query by bin range
        Index("ix_sampleID", sampleID),
        {
            "mysql_engine": "MyISAM",
            "mysql_charset": "utf8"
        }
    )

    @classmethod
    def is_unique(cls, session, regionID,
        sampleID, experimentID, sourceID,
        histone_type):

        q = session.query(cls).\
            filter(
                cls.regionID == regionID,
                cls.sampleID == sampleID,
                cls.experimentID == experimentID,
                cls.sourceID == sourceID,
                cls.histone_type == histone_type
            )

        return len(q.all()) == 0

    @classmethod
    def select_unique(cls, session, regionID,
        sampleID, experimentID, sourceID,
        histone_type):

        q = session.query(cls).\
            filter(
                cls.regionID == regionID,
                cls.sampleID == sampleID,
                cls.experimentID == experimentID,
                cls.sourceID == sourceID,
                cls.histone_type == histone_type
            )

        return q.first()

    def __repr__(self):

        return "<%s(%s, %s, %s, %s, %s, %s)>" % \
            (
                self.__tablename__,
                "uid={}".format(self.uid),
                "regionID={}".format(self.regionID),
                "sampleID={}".format(self.sampleID),
                "experimentID={}".format(self.experimentID),
                "sourceID={}".format(self.sourceID),
                "histone_type={}".format(self.histone_type)
            )