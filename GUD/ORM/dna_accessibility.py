from binning import (containing_bins, contained_bins)
from sqlalchemy import (UniqueConstraint, Index)
from sqlalchemy.dialects import mysql
from .base import Base
from .experiment import Experiment
from .region import Region
from .sample import Sample
from .source import Source
from .genomic_feature import GenomicFeature
from .genomicFeatureMixin2 import GFMixin2
from sqlalchemy.ext.declarative import declared_attr

class DNAAccessibility(GFMixin2, Base):

    __tablename__ = "dna_accessibility"

    @declared_attr
    def __table_args__(cls):
        return (
            UniqueConstraint(
                cls.region_id,
                cls.sample_id,
                cls.experiment_id,
                cls.sample_id
            ),
            Index("ix_regionID", cls.region_id),  # query by bin range
            Index("ix_sourceID", cls.sample_id),
            Index("ix_experimentID", cls.experiment_id),
            Index("ix_sampleID", cls.sample_id),
            {
                "mysql_engine": "InnoDB",
                "mysql_charset": "utf8"
            }
        )

    @classmethod
    def select_unique(cls, session, regionID,
                      sampleID, experimentID, sourceID):

        q = session.query(cls).\
            filter(
                cls.regionID == regionID,
                cls.sampleID == sampleID,
                cls.experimentID == experimentID,
                cls.sourceID == sourceID
        )

        return q.first()

    @classmethod
    def as_genomic_feature(self, feat):
        # Define qualifiers
        qualifiers = {
            "uid": feat.DNAAccessibility.uid,
            "source": feat.Source.name,
            "sample": feat.Sample.name,
            "experiment": feat.Experiment.name,
        }
        return GenomicFeature(
            feat.Region.chrom,
            int(feat.Region.start),
            int(feat.Region.end),
            strand=feat.Region.strand,
            feat_id="%s_%s" % (self.__tablename__, feat.DNAAccessibility.uid),
            qualifiers=qualifiers)
