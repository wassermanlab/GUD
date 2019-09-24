from binning import (
    containing_bins,
    contained_bins
)
from sqlalchemy import (
    Column,
    Index,
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


class TFBinding(GFMixin2, Base):

    __tablename__ = "tf_binding"

    tf = Column("tf", String(25), nullable=False)

    @declared_attr
    def __table_args__(cls):
        return (
            UniqueConstraint(
                cls.region_id,
                cls.sample_id,
                cls.experiment_id,
                cls.sample_id,
                cls.tf
            ),
            Index("ix_regionID", cls.region_id),  # query by bin range
            Index("ix_sampleID", cls.sample_id),
            {
                "mysql_engine": "InnoDB",
                "mysql_charset": "utf8"
            }
        )

    @classmethod
    def select_by_tf(cls, query, tf):
        """
        Query objects by sources.
        """
        q = query.filter(cls.tf.in_(tf))

        return q

    @classmethod
    def is_unique(cls, session, regionID,
                  sampleID, experimentID, sourceID,
                  tf):

        q = session.query(cls).\
            filter(cls.regionID == regionID, cls.sampleID == sampleID,
                   cls.experimentID == experimentID, cls.sourceID == sourceID,
                   cls.tf == tf)

        return len(q.all()) == 0

    @classmethod
    def select_unique(cls, session, regionID,
                      sampleID, experimentID, sourceID,
                      tf):

        q = session.query(cls).\
            filter(cls.regionID == regionID, cls.sampleID == sampleID,
                   cls.experimentID == experimentID, cls.sourceID == sourceID,
                   cls.tf == tf)

        return q.first()

    @classmethod
    def as_genomic_feature(self, feat):
        # Define qualifiers
        qualifiers = {
            "uid": feat.TFBinding.uid,
            "source": feat.Source.name,
            "sample": feat.Sample.name,
            "experiment": feat.Experiment.name,
            "tf": feat.TFBinding.tf
        }
        return GenomicFeature(
            feat.Region.chrom,
            int(feat.Region.start),
            int(feat.Region.end),
            strand=feat.Region.strand,
            feat_id="%s_%s" % (self.__tablename__, feat.TFBinding.uid),
            qualifiers=qualifiers)
