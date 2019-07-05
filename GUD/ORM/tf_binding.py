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
                "mysql_engine": "MyISAM",
                "mysql_charset": "utf8"
            }
        )

    @classmethod
    def select_by_tf(cls, session, chrom, start, end, tf, limit, offset):
        """
        Query objects by sources.
        """
        bins = Region._compute_bins(start, end)

        q = session.query(cls, Region, Source, Sample, Experiment)\
            .filter(Region.uid == cls.region_id, Source.uid == cls.source_id,
                    Sample.uid == cls.sample_id, Experiment.uid == cls.experiment_id)\
            .filter(Region.chrom == chrom,
                    Region.start < end,
                    Region.end > start)\
            .filter(Region.bin.in_(bins))
        
        res = []
        for i in q.all():
            if i.TFBinding.tf in tf:
                res.append(i)
        return (len(res), res[offset:offset+limit])

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
