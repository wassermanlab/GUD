from sqlalchemy import (Column, UniqueConstraint, Index)
from sqlalchemy.dialects import mysql
from .base import Base
from .experiment import Experiment
from .region import Region
from .sample import Sample
from .source import Source
from .genomicFeatureMixin2 import GFMixin2
from sqlalchemy.ext.declarative import declared_attr


class DNAAccessibility(GFMixin2, Base):
    # table declerations
    __tablename__ = "dna_accessibility"
    score = Column("score", mysql.FLOAT)
    peak = Column("peak", mysql.INTEGER)

    @declared_attr
    def __table_args__(cls):
        return (
        UniqueConstraint(cls.region_id, cls.sample_id, cls.experiment_id,
                         cls.source_id, cls.peak),
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
            "source": feat.sourceName,
            "sample": feat.sampleName,
            "experiment": feat.experimentName,
            "score": feat.DNAAccessibility.score,
            "peak": feat.DNAAccessibility.peak
        }
        genomic_feature = super().as_genomic_feature(feat)
        genomic_feature.qualifiers = qualifiers
        return genomic_feature
