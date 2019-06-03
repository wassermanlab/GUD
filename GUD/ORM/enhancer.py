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
from .genomic_feature import GenomicFeature
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


    @classmethod
    def select_by_sample(cls, session,
        sample, as_genomic_feature=False):
        """
        Query objects by sample.
        """

        q = session.query(
            cls,
            Region,
            Sample,
            Experiment,
            Source,
            
        )\
            .join()\
            .filter(
                Region.uid == cls.regionID,
                Sample.uid == cls.sampleID,
                Experiment.uid == cls.experimentID,
                Source.uid == cls.sourceID
            )\
            .filter(
                Sample.name == sample
            )

        if as_genomic_feature:

            feats = []

            # For each feature...
            for feat in q.all():
                feats.append(
                    cls.__as_genomic_feature(feat)
                )

            return feats
    
        return q.all()

    def __as_genomic_feature(feat):

        # Define qualifiers
        qualifiers = {
            "experiment": feat.Experiment.name,
            "sample": feat.Sample.name,
            "source" : feat.Source.name,            
        }

        return GenomicFeature(
            feat.Region.chrom,
            int(feat.Region.start),
            int(feat.Region.end),
            strand = feat.Region.strand,
            feat_type = "Enhancer",
            feat_id = "%s|%s" % (
                qualifiers["source"],
                qualifiers["sample"].replace(" ", "_")
            ),
            qualifiers = qualifiers
        )

    def __repr__(self):

        return "<Enhancer(%s, %s, %s, %s, %s)>" % \
            (
                "uid={}".format(self.uid),
                "regionID={}".format(self.regionID),
                "sampleID={}".format(self.sampleID),
                "experimentID={}".format(self.experimentID),
                "sourceID={}".format(self.sourceID)
            )