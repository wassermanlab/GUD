from sqlalchemy import (Column, Index, String, UniqueConstraint)
from sqlalchemy.dialects import mysql
from .base import Base
from .experiment import Experiment
from .region import Region
from .sample import Sample
from .source import Source
from .genomicFeatureMixin2 import GFMixin2
from sqlalchemy.ext.declarative import declared_attr


class TAD(GFMixin2, Base):
    # table declerations
    __tablename__ = "tads"
    hierarchical_level = Column("hierarchical_level", String(25))

    @declared_attr
    def __table_args__(cls):
        return (
            UniqueConstraint(cls.region_id, cls.sample_id,
                             cls.experiment_id, cls.sample_id),
            Index("ix_join", cls.region_id, cls.sample_id,
                  cls.experiment_id, cls.source_id),
            {"mysql_engine": "InnoDB", "mysql_charset": "utf8"}
        )

    # class methods
    @classmethod
    def is_unique(cls, session, regionID, sampleID, experimentID,
                  sourceID):

        q = session.query(cls).\
            filter(cls.region_id == regionID, cls.sample_id == sampleID,
                   cls.experiment_id == experimentID, cls.source_id == sourceID)

        return len(q.all()) == 0

    def as_genomic_feature(self, feat):

        # Define qualifiers
        qualifiers = {
            "uid": feat.TAD.uid,
            "source": feat.Source.name,
            "sample": feat.Sample.name,
            "experiment": feat.Experiment.name,
            "hierarchical_level": feat.TAD.hierarchical_level
        }
        genomic_feature = super().as_genomic_feature(feat)
        genomic_feature.qualifiers = qualifiers
        return genomic_feature
