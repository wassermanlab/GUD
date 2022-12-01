from sqlalchemy import (Column, ForeignKey)
from sqlalchemy.dialects import mysql
from .region import Region
from .source import Source
from .sample import Sample
from .experiment import Experiment
from sqlalchemy.ext.declarative import declared_attr
from .genomic_feature import GenomicFeature
from .genomicFeatureMixin1 import GFMixin1

class GFMixin2(GFMixin1):
    @declared_attr
    def sample_id(cls):
        return Column("sampleID", ForeignKey("samples.uid"),
                      nullable=False)

    @declared_attr
    def experiment_id(cls):
        return Column("experimentID", ForeignKey("experiments.uid"),
                      nullable=False)

    # methods###################
    @classmethod
    def get_unique_sample_names(cls, session):
        sample_ids = session.query(cls.sample_id).distinct().all()
        sample_ids = [s[0] for s in sample_ids]
        sample_names = session.query(Sample.name).filter(Sample.uid.in_(sample_ids)).distinct().all()
        sample_names = [s[0] for s in sample_names]
        return sample_names

    @classmethod
    def get_unique_experiment_names(cls, session):
        experiment_ids = session.query(cls.experiment_id).distinct().all()
        experiment_ids = [s[0] for s in experiment_ids]
        experiment_names = session.query(Experiment.name).filter(Experiment.uid.in_(experiment_ids)).distinct().all()
        experiment_names = [s[0] for s in experiment_names]
        return experiment_names
    
    @classmethod
    def make_query(cls, session, query):
        if (query is not None):
            return query
            # Region, Source, Sample, Experiment
        q = session.query(cls, Region.chrom, Region.start, Region.end, Source.name.label('sourceName'), Sample.name.label('sampleName'), Experiment.name.label('experimentName'))\
            .prefix_with("STRAIGHT_JOIN")\
            .join(Region, Region.uid == cls.region_id)\
            .join(Source, Source.uid == cls.source_id)\
            .join(Experiment, Experiment.uid == cls.experiment_id)\
            .join(Sample, Sample.uid == cls.sample_id)
        return q
        
    @classmethod
    def select_by_samples(cls, session, query, samples):
        """
        filter query by samples.
        """
        sample_uids = session.query(Sample).filter(Sample.name.in_(samples)).all()
        sample_uids = [t.uid for t in sample_uids]
        q = cls.make_query(session, query)
        q = q.filter(Sample.uid.in_(sample_uids))
        return q

    @classmethod
    def select_by_experiments(cls, session, query, experiments):
        """
        filter query by experiments.
        """
        experiment_uids = session.query(Experiment).filter(Experiment.name.in_(experiments)).all()
        experiment_uids = [t.uid for t in experiment_uids]
        q = cls.make_query(session, query)
        q = q.filter(Experiment.uid.in_(experiment_uids))
        return q

    @classmethod
    def is_unique(cls, session, regionID, sampleID, experimentID, sourceID):
        """
        Check by unique condition. Minimum condition is is regionID and sourceID.
        Additional conditions require overriding of method in child class.
        """
        q = session.query(cls).filter(cls.region_id == regionID,
                                      cls.sample_id == sampleID,
                                      cls.experiment_id == experimentID,
                                      cls.source_id == sourceID)
        q = q.all()
        return len(q) == 0