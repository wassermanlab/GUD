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
from .genomic_feature import GenomicFeature
from .region import Region
from .sample import Sample
from .source import Source


class TAD(Base):

    __tablename__ = "tads"

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

    restriction_enzyme = Column(
        "restriction_enzyme",
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
            restriction_enzyme

        ),
        Index("ix_regionID", regionID),  # query by bin range
        Index("ix_sampleID", sampleID),
        {
            "mysql_engine": "MyISAM",
            "mysql_charset": "utf8"
        }
    )

    @classmethod
    def is_unique(cls, session, regionID,
                  sampleID, experimentID, sourceID,
                  restriction_enzyme):

        q = session.query(cls).\
            filter(
                cls.regionID == regionID,
                cls.sampleID == sampleID,
                cls.experimentID == experimentID,
                cls.sourceID == sourceID,
                cls.restriction_enzyme == restriction_enzyme
        )

        return len(q.all()) == 0

    @classmethod
    def select_unique(cls, session, regionID,
                      sampleID, experimentID, sourceID,
                      restriction_enzyme):

        q = session.query(cls).\
            filter(
                cls.regionID == regionID,
                cls.sampleID == sampleID,
                cls.experimentID == experimentID,
                cls.sourceID == sourceID,
                cls.restriction_enzyme == restriction_enzyme
        )

        return q.first()

    @classmethod
    def select_by_location(cls, session,
                           chrom, start, end, samples=[],
                           as_genomic_feature=False):
        """
        Query objects by genomic location.
        """

        bins = Region._compute_bins(start, end)

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
                Region.chrom == chrom,
                Region.start < end,
                Region.end > start
        )\
            .filter(Region.bin.in_(bins))

        if samples:
            q = q.filter(Sample.name.in_(samples))

        if as_genomic_feature:

            feats = []

            # For each feature...
            for feat in q.all():
                feats.append(
                    cls.__as_genomic_feature(feat)
                )

            return feats

        return q.all()

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

    def __as_genomic_feature(self, feat):

        # Define qualifiers
        qualifiers = {
            "uid": feat.TAD.uid,
            "regionID": feat.TAD.regionID,
            "sampleID": feat.TAD.sampleID,
            "experimentID": feat.TAD.experimentID,
            "sourceID": feat.TAD.sourceID,
            "restriction_enzyme": feat.TAD.restriction_enzyme,
            "experiment": feat.Experiment.name,
            "sample": feat.Sample.name,
            "source": feat.Source.name,
        }

        return GenomicFeature(
            feat.Region.chrom,
            int(feat.Region.start),
            int(feat.Region.end),
            strand=feat.Region.strand,
            feat_type="TAD",
            feat_id="%s|%s|%s" % (
                qualifiers["source"],
                qualifiers["restriction_enzyme"],
                qualifiers["sample"].replace(" ", "_")
            ),
            qualifiers=qualifiers
        )

    def __repr__(self):

        return "<TAD(%s, %s, %s, %s, %s, %s)>" % \
            (
                "uid={}".format(self.uid),
                "regionID={}".format(self.regionID),
                "sampleID={}".format(self.sampleID),
                "experimentID={}".format(self.experimentID),
                "sourceID={}".format(self.sourceID),
                "restriction_enzyme={}".format(self.restriction_enzyme)
            )
