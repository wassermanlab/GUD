from sqlalchemy import (
    Column,
    Index,
    Integer,
    String,
    ForeignKey,
    UniqueConstraint
) 
import re 
from sqlalchemy.dialects import mysql
from .base import Base
from .experiment import Experiment
from .gene import Gene
from .region import Region
from .source import Source
from .genomic_feature import GenomicFeature
from .genomicFeatureMixin2 import GFMixin2
from sqlalchemy.ext.declarative import declared_attr


class TSS(GFMixin2, Base):

    __tablename__ = "transcription_start_sites"

    gene = Column("gene", String(75), ForeignKey("genes.name2"))

    tss = Column("tss", mysql.INTEGER(unsigned=True))

    @declared_attr
    def sample_id(cls):
        return Column("sampleIDs", mysql.LONGBLOB, nullable=False)

    avg_expression_levels = Column("avg_expression_levels", mysql.LONGBLOB,
                                   nullable=False)

    @declared_attr
    def __table_args__(cls):
        return (
            # multiple TSSs might overlap:
            # e.g. p16@IGF2,p1@INS-IGF2,p1@INS
            UniqueConstraint(
                cls.region_id,
                cls.sample_id,
                cls.experiment_id,
                cls.gene,
                cls.tss
            ),
            Index("ix_regionID", cls.region_id),  # query by bin range
            Index("ix_gene_tss", cls.gene, cls.tss),
            {
                "mysql_engine": "InnoDB",
                "mysql_charset": "utf8"
            }
        )

    @classmethod
    def is_unique(cls, session, regionID, sourceID,
                  experimentID, gene, tss):

        q = session.query(cls)\
            .filter(cls.regionID == regionID, cls.sourceID == sourceID,
                    cls.experimentID == experimentID, cls.gene == gene, cls.tss == tss)

        return len(q.all()) == 0

    @classmethod
    def select_unique(cls, session, regionID, sourceID, experimentID, gene, tss):

        q = session.query(cls)\
            .filter(cls.regionID == regionID, cls.sourceID == sourceID,
                    cls.experimentID == experimentID, cls.gene == gene,
                    cls.tss == tss)

        return q.first()

    @classmethod
    def select_by_overlapping_location(cls, session, chrom, start, end):
        """
        Query objects by genomic location.
        """
        bins = Region._compute_bins(start, end)

        q = session.query(cls, Region, Source, Experiment)\
            .filter(Region.uid == cls.region_id, Source.uid == cls.source_id,
                    Experiment.uid == cls.experiment_id)\
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

        q = session.query(cls, Region, Source, Experiment)\
            .filter(Region.uid == cls.region_id, Source.uid == cls.source_id,
                    Experiment.uid == cls.experiment_id)\
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

        q = session.query(cls, Region, Source, Experiment)\
            .filter(Region.uid == cls.region_id, Source.uid == cls.source_id,
                    Experiment.uid == cls.experiment_id)\
            .filter(Region.bin.in_(bins))\
            .filter(Region.uid == cls.region_id)\
            .filter(Region.chrom == chrom,
                    Region.start == start,
                    Region.end == end)
        return q

    @classmethod
    def select_by_uids(cls, session, query, uids):
        """
        Query objects by uids.
        """
        q = query.filter(cls.uid.in_(uids))

        return q

    @classmethod
    def select_by_genes(cls, query, genes):
        """
        select by genes
        """

        q = query.filter(cls.gene.in_(genes))

        return q

    @classmethod  # TODO
    def select_all_genic_tss(cls, session):
        """
        Query all objects associated with a gene.
        """
        q = session.query(cls, Region, Source, Experiment).\
            filter(Region.uid == cls.region_id, Source.uid == cls.source_id,
                   Experiment.uid == cls.experiment_id).filter(cls.gene != None)

        return q

    @classmethod
    def select_by_sources(cls, session, query, sources):
        """
        Query objects by sources.
        """
        s = [value for value, in session.query(
            Source.uid).filter(Source.name.in_(sources)).all()]
        q = query.filter(Source.uid.in_(s))

        return q

    @classmethod
    def select_by_samples(cls, session, query, samples):
        """
        Query objects by sample name.
        """
        return query

    @classmethod  # TODO
    def select_by_experiments(cls, session, query, experiments):
        """
        Query objects by experiment name.
        """

        s = [value for value, in session.query(Experiment.uid).filter(
            Experiment.name.in_(experiments)).all()]
        q = query.filter(Experiment.uid.in_(s))

        return q

    @classmethod
    def as_genomic_feature(self, feat):

        # Initialize
        isfloat = re.compile("\d+(\.\d+)?")
        sampleIDs = []
        avg_expression_levels = []

        # For each exon start...
        for i in str(feat.TSS.sample_id).split(","):
            if i.isdigit():
                sampleIDs.append(int(i))

        # For each exon end...
        for i in str(feat.TSS.avg_expression_levels).split(","):
            if isfloat.match(i):
                avg_expression_levels.append(float(i))

        # Define qualifiers
        qualifiers = {
            "uid": feat.TSS.uid,
            "gene": feat.TSS.gene,
            "tss": feat.TSS.tss,
            "sampleIDs": feat.TSS.sample_id,
            "avg_expression_levels": feat.TSS.avg_expression_levels,
            "experiment": feat.Experiment.name,
            "source": feat.Source.name,
        }

        return GenomicFeature(
            feat.Region.chrom,
            int(feat.Region.start),
            int(feat.Region.end),
            strand=feat.Region.strand,
            feat_type="TSS",
            feat_id="%s_%s" % (self.__tablename__, feat.TSS.uid),
            qualifiers=qualifiers
        )
