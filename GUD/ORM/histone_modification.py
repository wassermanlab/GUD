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


class HistoneModification(GFMixin2, Base):

    __tablename__ = "histone_modifications"

    histone_type = Column("histone_type", String(25), nullable=False)

    
    
    @declared_attr
    def __table_args__(cls):
        return (
            UniqueConstraint(
                cls.region_id,
                cls.sample_id,
                cls.experiment_id,
                cls.sample_id,
                cls.histone_type

            ),
            Index("ix_regionID", cls.region_id),  # query by bin range
            Index("ix_sampleID", cls.sample_id),
            {
                "mysql_engine": "InnoDB",
                "mysql_charset": "utf8"
            }
        )

    @classmethod
    def select_by_histone_type(cls, query, histone_type):
        """
        Query objects by sources.
        """

        q = query.filter(cls.histone_type.in_(histone_type))

        return q

    @classmethod
    def is_unique(cls, session, regionID, sampleID, experimentID, sourceID,
                  histone_type):

        q = session.query(cls).\
            filter(cls.regionID == regionID, cls.sampleID == sampleID,
                   cls.experimentID == experimentID, cls.sourceID == sourceID,
                   cls.histone_type == histone_type)

        return len(q.all()) == 0

    @classmethod
    def select_unique(cls, session, regionID, sampleID, experimentID, sourceID,
                      histone_type):

        q = session.query(cls).\
            filter(cls.regionID == regionID, cls.sampleID == sampleID,
                   cls.experimentID == experimentID, cls.sourceID == sourceID,
                   cls.histone_type == histone_type)

        return q.first()

    @classmethod
    def as_genomic_feature(self, feat):
        # Define qualifiers
        qualifiers = {
            "uid": feat.HistoneModification.uid,
            "histone_type": feat.HistoneModification.histone_type,
            "source": feat.Source.name,
            "sample": feat.Sample.name,
            "experiment": feat.Experiment.name,
        }
        return GenomicFeature(
            feat.Region.chrom,
            int(feat.Region.start),
            int(feat.Region.end),
            strand = feat.Region.strand,
            feat_id = "%s_%s"%(self.__tablename__, feat.HistoneModification.uid),
            qualifiers = qualifiers)

