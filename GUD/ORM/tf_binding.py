from sqlalchemy import (Column, Index, String, UniqueConstraint)
from sqlalchemy.dialects import mysql
from .base import Base
from .experiment import Experiment
from .region import Region
from .sample import Sample
from .source import Source
from .genomicFeatureMixin2 import GFMixin2
from sqlalchemy.ext.declarative import declared_attr


class TFBinding(GFMixin2, Base):
    # table declerations 
    __tablename__ = "tf_binding"
    tf = Column("tf", String(25), nullable=False)
    peak =  Column("peak", mysql.LONGBLOB, nullable=False)

    @declared_attr
    def __table_args__(cls):
        return (
        UniqueConstraint(cls.region_id, cls.sample_id,
                         cls.experiment_id, cls.source_id, cls.tf),
        Index("ix_join", cls.region_id, cls.sample_id, cls.experiment_id, cls.source_id),
        {"mysql_engine": "InnoDB", "mysql_charset": "utf8"}
    )

    # class methods 
    @classmethod
    def select_by_tf(cls, session, query, tf):
        """
        Query objects by sources.
        """
        q = cls.make_query(session, query)
        q = q.filter(cls.tf.in_(tf))

        return q

    @classmethod
    def is_unique(cls, session, regionID,
                  sampleID, experimentID, sourceID,
                  tf):

        q = session.query(cls).\
            filter(cls.region_id == regionID, cls.sample_id == sampleID,
                   cls.experiment_id == experimentID, cls.source_id == sourceID,
                   cls.tf == tf)

        return len(q.all()) == 0

    @classmethod
    def as_genomic_feature(self, feat):
        # Define qualifiers
        qualifiers = {
            "uid": feat.TFBinding.uid,
            "source": feat.Source.name,
            "sample": feat.Sample.name,
            "experiment": feat.Experiment.name,
            "tf": feat.TFBinding.tf,
            "peak": feat.HistoneModification.peak
        }
        genomic_feature = super().as_genomic_feature(feat)
        genomic_feature.qualifiers = qualifiers
        return genomic_feature
