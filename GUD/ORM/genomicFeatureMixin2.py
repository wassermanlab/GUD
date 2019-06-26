from sqlalchemy import (
    Column,
    ForeignKey,
    Integer,
    PrimaryKeyConstraint
)

from sqlalchemy.dialects import mysql
from .base import Base
from .genomic_feature import GenomicFeature
from .region import Region
from .source import Source
from .sample import Sample
from .experiment import Experiment
from sqlalchemy.ext.declarative import declared_attr
from .genomicFeatureMixin1 import GFMixin1


class GFMixin2(GFMixin1):

    @declared_attr
    def sample_id(cls):
        return Column("sampleID", Integer, ForeignKey("samples.uid"),
                      nullable=False)

    @declared_attr
    def experiment_id(cls):
        return Column("experimentID", Integer, ForeignKey("experiments.uid"),
                      nullable=False)

    @classmethod
    def select_by_uids(cls, session, uids, limit, offset):
        """
        Query objects by uids.
        """
        q = session.query(cls, Region, Source, Sample, Experiment).\
            join()\
            .filter(Region.uid == cls.region_id, Source.uid == cls.source_id,
                    Sample.uid == cls.sample_id, Experiment.uid == cls.experiment_id)\
            .filter(cls.uid.in_(uids))

        return (q.count(), q.offset(offset).limit(limit))

    @classmethod
    def select_by_sources(cls, session, sources, limit, offset):
        """
        Query objects by sources.
        """
        q = session.query(cls, Region, Source, Sample, Experiment).\
            join()\
            .filter(Region.uid == cls.region_id, Source.uid == cls.source_id,
                    Sample.uid == cls.sample_id, Experiment.uid == cls.experiment_id)\
            .filter(Source.name.in_(sources))

        return (q.count(), q.offset(offset).limit(limit))

    @classmethod
    def select_by_samples(cls, session, samples, limit, offset):
        """
        Query objects by sample name.
        """
        q = session.query(cls, Region, Source, Sample, Experiment).\
            join()\
            .filter(Region.uid == cls.region_id, Source.uid == cls.source_id,
                    Sample.uid == cls.sample_id, Experiment.uid == cls.experiment_id)\
            .filter(Sample.name.in_(samples))

        return (q.count(), q.offset(offset).limit(limit))

    @classmethod
    def select_by_experiments(cls, session, experiments, limit, offset):
        """
        Query objects by experiment name.
        """
        q = session.query(cls, Region, Source, Sample, Experiment).\
            join()\
            .filter(Region.uid == cls.region_id, Source.uid == cls.source_id,
                    Sample.uid == cls.sample_id, Experiment.uid == cls.experiment_id)\
            .filter(Experiment.name.in_(experiments))

        return (q.count(), q.offset(offset).limit(limit))

    @classmethod
    def select_by_location(cls, session, chrom, start, end, limit, offset):
        """
        Query objects by genomic location.
        """
        bins = Region._compute_bins(start, end)

        q = session.query(cls, Region, Source, Sample, Experiment)\
            .join()\
            .filter(Region.uid == cls.region_id, Source.uid == cls.source_id,
                    Sample.uid == cls.sample_id, Experiment.uid == cls.experiment_id)\
            .filter(Region.chrom == chrom,
                    Region.start < end,
                    Region.end > start)\
            .filter(Region.bin.in_(bins))

        return (q.count(), q.offset(offset).limit(limit))

    @classmethod
    def select_by_exact_location(cls, session, chrom, start, end, limit, offset):
        """
        Query objects by genomic location.
        """
        bins = Region._compute_bins(start, end)

        q = session.query(cls, Region, Source, Sample, Experiment)\
            .join()\
            .filter(Region.bin.in_(bins))\
            .filter(Region.uid == cls.region_id, Source.uid == cls.source_id,
                    Sample.uid == cls.sample_id, Experiment.uid == cls.experiment_id)\
            .filter(Region.chrom == chrom,
                    Region.start == start,
                    Region.end == end)

        return (q.count(), q.offset(offset).limit(limit))

    @classmethod
    def is_unique(cls, session, regionID,
                  sampleID, experimentID, sourceID):

        q = session.query(cls).\
            filter(cls.regionID == regionID, cls.sampleID == sampleID,
                   cls.experimentID == experimentID, cls.sourceID == sourceID)

        return len(q.all()) == 0

    @classmethod
    def select_unique(cls, session, regionID,
                      sampleID, experimentID, sourceID):

        q = session.query(cls).\
            filter(cls.regionID == regionID, cls.sampleID == sampleID,
                   cls.experimentID == experimentID, cls.sourceID == sourceID)

        return q.first()
