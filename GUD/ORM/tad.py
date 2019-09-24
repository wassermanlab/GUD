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
from .genomic_feature import GenomicFeature
from .genomicFeatureMixin2 import GFMixin2
from sqlalchemy.ext.declarative import declared_attr


class TAD(GFMixin2, Base):

    __tablename__ = "tads"

    restriction_enzyme = Column(
        "restriction_enzyme", String(25), nullable=False)

    @declared_attr
    def __table_args__(cls):
        return (
            UniqueConstraint(
                cls.region_id,
                cls.sample_id,
                cls.experiment_id,
                cls.sample_id,
                cls.restriction_enzyme

            ),
            Index("ix_regionID", cls.region_id),  # query by bin range
            Index("ix_sampleID", cls.sample_id),
            {
                "mysql_engine": "InnoDB",
                "mysql_charset": "utf8"
            }
        )

    @classmethod
    def select_by_restriction_enzymes(cls, query, restriction_enzymes):
        """
        Query objects by sources.
        """

        q = query.filter(cls.restriction_enzyme.in_(restriction_enzymes))

        return q

    @classmethod
    def is_unique(cls, session, regionID,
                  sampleID, experimentID, sourceID,
                  restriction_enzyme):

        q = session.query(cls).\
            filter(cls.regionID == regionID, cls.sampleID == sampleID,
                   cls.experimentID == experimentID, cls.sourceID == sourceID,
                   cls.restriction_enzyme == restriction_enzyme)

        return len(q.all()) == 0

    @classmethod
    def select_unique(cls, session, regionID, sampleID, experimentID, sourceID,
                      restriction_enzyme):

        q = session.query(cls).\
            filter(cls.regionID == regionID, cls.sampleID == sampleID,
                   cls.experimentID == experimentID, cls.sourceID == sourceID,
                   cls.restriction_enzyme == restriction_enzyme)

        return q.first()

    def as_genomic_feature(self, feat):

        # Define qualifiers
        qualifiers = {
            "uid": feat.TAD.uid,
            "source": feat.Source.name,
            "sample": feat.Sample.name,
            "experiment": feat.Experiment.name,
            "restriction_enzyme": feat.TAD.restriction_enzyme,
        }

        return GenomicFeature(
            feat.Region.chrom,
            int(feat.Region.start),
            int(feat.Region.end),
            strand=feat.Region.strand,
            feat_type="TAD",
            feat_id="%s_%s" % (self.__tablename__, feat.TAD.uid),
            qualifiers=qualifiers
        )
