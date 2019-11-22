from sqlalchemy import (Index, UniqueConstraint)
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
        UniqueConstraint(cls.region_id, cls.sample_id,
                         cls.experiment_id, cls.source_id),
        Index("ix_join", cls.region_id, cls.sample_id,
              cls.experiment_id, cls.source_id),
        {"mysql_engine": "InnoDB", "mysql_charset": "utf8"}
    )

    # class methods
    @classmethod
    def as_genomic_feature(self, feat):
        # Define qualifiers
        qualifiers = {
            "uid": feat.DNAAccessibility.uid,
            "source": feat.Source.name,
            "sample": feat.Sample.name,
            "experiment": feat.Experiment.name,
        }
        genomic_feature = super().as_genomic_feature(feat)
        genomic_feature.qualifiers = qualifiers
        return genomic_feature
