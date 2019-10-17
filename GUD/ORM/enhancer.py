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
from .genomic_feature import GenomicFeature
from .genomicFeatureMixin2 import GFMixin2
from sqlalchemy.ext.declarative import declared_attr


class Enhancer(GFMixin2, Base):

    __tablename__ = "enhancers"

    @declared_attr
    def __table_args__(cls):
        return (
        UniqueConstraint(GFMixin2.region_id, GFMixin2.sample_id,
                         GFMixin2.experiment_id, GFMixin2.source_id),
        Index("ix_join", GFMixin2.region_id, GFMixin2.sample_id,
              GFMixin2.experiment_id, GFMixin2.source_id),
        {"mysql_engine": "InnoDB", "mysql_charset": "utf8"}
    )

    def as_genomic_feature(self, feat):
        qualifiers = {
            "uid": feat.Enhancer.uid,
            "source": feat.Source.name,
            "sample": feat.Sample.name,
            "experiment": feat.Experiment.name,
        }

        return GenomicFeature(
            feat.Region.chrom,
            int(feat.Region.start),
            int(feat.Region.end),
            strand=feat.Region.strand,
            feat_type=self.__tablename__,
            feat_id="%s_%s" %
            (
                self.__tablename__,
                feat.Enhancer.uid
            ),
            qualifiers=qualifiers
        )
