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
import sys


class GFMixin2(GFMixin1):

    @declared_attr
    def sample_id(cls):
        return Column("sampleID", mysql.INTEGER(unsigned=True), ForeignKey("samples.uid"),
                      nullable=False)

    @declared_attr
    def experiment_id(cls):
        return Column("experimentID", mysql.INTEGER(unsigned=True), ForeignKey("experiments.uid"),
                      nullable=False)

    @classmethod
    def select_by_overlapping_location(cls, session, chrom, start, end):
        """
        Query objects by genomic location.
        """
        bins = Region._compute_bins(start, end)

        q = session.query(cls, Region, Source, Sample, Experiment)\
            .filter(Region.uid == cls.region_id, Source.uid == cls.source_id,
                    Sample.uid == cls.sample_id, Experiment.uid == cls.experiment_id)\
            .filter(Region.chrom == chrom,
                    Region.start < end,
                    Region.end > start)\
            .filter(Region.bin.in_(bins))

        return q
    
    @classmethod
    def select_by_within_location(cls, session, chrom, start, end):
        """
        Query objects by genomic location.
        """
        bins = Region._compute_bins(start, end)

        q = session.query(cls, Region, Source, Sample, Experiment)\
            .filter(Region.uid == cls.region_id, Source.uid == cls.source_id,
                    Sample.uid == cls.sample_id, Experiment.uid == cls.experiment_id)\
            .filter(Region.chrom == chrom,
                    Region.start > start,
                    Region.end < end)\
            .filter(Region.bin.in_(bins))
        
        return q

    @classmethod
    def select_by_exact_location(cls, session, chrom, start, end):
        """
        Query objects by genomic location.
        """
        bins = Region._compute_bins(start, end)

        q = session.query(cls, Region, Source, Sample, Experiment)\
            .filter(Region.uid == cls.region_id, Source.uid == cls.source_id,
                    Sample.uid == cls.sample_id, Experiment.uid == cls.experiment_id)\
            .filter(Region.bin.in_(bins))\
            .filter(Region.uid == cls.region_id)\
            .filter(Region.chrom == chrom,
                    Region.start == start,
                    Region.end == end)
        return q

    @classmethod
    def select_by_location(cls, session, chrom, start, end, location):
        
        if location == 'exact':
            q = cls.select_by_exact_location(session, chrom, start, end)
        elif location == 'within':
            q = cls.select_by_within_location(session, chrom, start, end)
        elif location == 'overlapping':
            q = cls.select_by_overlapping_location(session, chrom, start, end)
        else: 
            return False
        
        return q

    @classmethod
    def select_by_uids(cls, session, query, uids):
        """
        Query objects by uids.
        """
        q = query.filter(cls.uid.in_(uids))

        return q

    @classmethod
    def select_by_sources(cls, session, query, sources):
        """
        Query objects by sources.
        """
        s = [value for value, in session.query(Source.uid).filter(Source.name.in_(sources)).all()]
        q = query.filter(Source.uid.in_(s))

        return q
            # .filter(Source.name.in_(sources))
        # return (q.count(), q.offset(offset).limit(limit))

    @classmethod
    def select_by_samples(cls, session, query, samples):
        """
        Query objects by sample name.
        """
        s = [value for value, in session.query(Sample.uid).filter(Sample.name.in_(samples)).all()]
        q = query.filter(Sample.uid.in_(s))

        return q
       
    @classmethod
    def select_by_experiments(cls, session, query, experiments):
        """
        Query objects by experiment name.
        """
        s = [value for value, in session.query(Experiment.uid).filter(Experiment.name.in_(experiments)).all()]
        q = query.filter(Experiment.uid.in_(s))
        
        return q

    @classmethod
    def is_unique(cls, session, regionID,
                  sampleID, experimentID, sourceID):

        q = session.query(cls).\
            filter(cls.region_id == regionID, cls.sample_id == sampleID,
                   cls.experiment_id == experimentID, cls.source_id == sourceID)

        return len(q.all()) == 0

    @classmethod
    def select_unique(cls, session, regionID,
                      sampleID, experimentID, sourceID):

        q = session.query(cls).\
            filter(cls.region_id == regionID, cls.sample_id == sampleID,
                   cls.experiment_id == experimentID, cls.source_id == sourceID)

        return q.first()
