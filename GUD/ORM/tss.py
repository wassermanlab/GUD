import re
from sqlalchemy import (
    and_,
    or_,
    Column,
    Index,
    Integer,
    PrimaryKeyConstraint,
    String,
    ForeignKey,
    UniqueConstraint
)
from sqlalchemy.dialects import mysql

from .base import Base
from .experiment import Experiment
from .gene import Gene
from .genomic_feature import GenomicFeature
from .region import Region
from .source import Source

class TSS(Base):

    __tablename__ = "transcription_start_sites"

    uid = Column(
        "uid",
        mysql.INTEGER(unsigned=True)
    )

    regionID = Column(
        "regionID",
        Integer,
        ForeignKey("regions.uid"),
        nullable=False
    )

    gene = Column(
        "gene",
        String(75),
        ForeignKey("genes.name2")
    )

    tss = Column(
        "tss",
        mysql.INTEGER(unsigned=True)
    )

    sampleIDs = Column(
        "sampleIDs",
        mysql.LONGBLOB,
        nullable=False
    )

    avg_expression_levels = Column(
        "avg_expression_levels",
        mysql.LONGBLOB, nullable=False
    )

    experimentID = Column(
        "experimentID",
        Integer,
        ForeignKey("experiments.uid"),
        nullable=False
    )

    sourceID = Column(
        "sourceID",
        Integer,
        ForeignKey("sources.uid"),
        nullable=False
    )

    __table_args__ = (
        PrimaryKeyConstraint(uid),
        # multiple TSSs might overlap:
        # e.g. p16@IGF2,p1@INS-IGF2,p1@INS
        UniqueConstraint(
            regionID,
            experimentID,
            sourceID,
            gene,
            tss
        ),
        Index("ix_regionID", regionID), # query by bin range
        Index("ix_gene_tss", gene, tss),
        {
            "mysql_engine": "MyISAM",
            "mysql_charset": "utf8"
        }
    )

    @classmethod
    def is_unique(cls, session, regionID, sourceID,
        experimentID, gene, tss):

        q = session.query(cls)\
            .filter(
                cls.regionID == regionID,
                cls.sourceID == sourceID,
                cls.experimentID == experimentID,
                cls.gene == gene,
                cls.tss == tss
            )

        return len(q.all()) == 0


    @classmethod
    def select_unique(cls, session, regionID,
        sourceID, experimentID, gene, tss):

        q = session.query(cls)\
            .filter(
                cls.regionID == regionID,
                cls.sourceID == sourceID,
                cls.experimentID == experimentID,
                cls.gene == gene,
                cls.tss == tss
            )

        return q.first()

    @classmethod
    def select_by_uid(cls, session, uid,
        as_genomic_feature=False):
        """
        Query objects by uid.
        """

        q = session.query(
                cls,
                Experiment,
                Region,
                Source
            )\
            .join()\
            .filter(
                Experiment.uid == cls.experimentID,
                Region.uid == cls.regionID,
                Source.uid == cls.sourceID
            ).filter(cls.uid == uid)

        if as_genomic_feature:
            return cls.__as_genomic_feature(
                q.first()
            )

        return q.first()

    @classmethod
    def select_by_uids(cls, session, uids=[],
        as_genomic_feature=False):
        """
        Query objects by multiple uids.
        If no uids are provided, return all
        objects.
        """

        q = session.query(
                cls,
                Experiment,
                Region,
                Source
            )\
            .join()\
            .filter(
                Experiment.uid == cls.experimentID,
                Region.uid == cls.regionID,
                Source.uid == cls.sourceID
            )

        if uids:
            q = q.filter(cls.uid.in_(uids))

        if as_genomic_feature:

            feats = []

            # For each feature...
            for feat in q.all():
                feats.append(
                    cls.__as_genomic_feature(feat)
                )

            return feats

        return q.all()

    @classmethod
    def select_by_gene(cls, session, gene,
        as_genomic_feature=False):
        """
        Query objects by uid.
        """

        q = session.query(
                cls,
                Experiment,
                Region,
                Source
            )\
            .join()\
            .filter(
                Experiment.uid == cls.experimentID,
                Region.uid == cls.regionID,
                Source.uid == cls.sourceID
            ).filter(cls.gene == gene)

        if as_genomic_feature:

            feats = []

            # For each feature...
            for feat in q.all():
                feats.append(
                    cls.__as_genomic_feature(feat)
                )

            return feats

        return q.all()

    @classmethod
    def select_by_genes(cls, session, genes=[],
        as_genomic_feature=False):
        """
        Query objects by multiple uids.
        If no uids are provided, return all
        objects.
        """

        q = session.query(
                cls,
                Experiment,
                Region,
                Source
            )\
            .join()\
            .filter(
                Experiment.uid == cls.experimentID,
                Region.uid == cls.regionID,
                Source.uid == cls.sourceID
            )

        if genes:
            q = q.filter(cls.gene.in_(genes))

        if as_genomic_feature:

            feats = []

            # For each feature...
            for feat in q.all():
                feats.append(
                    cls.__as_genomic_feature(feat)
                )

            return feats

        return q.all()

    @classmethod
    def select_all_genic_tss(cls, session,
        as_genomic_feature=False):
        """
        Query all objects associated with a gene.
        """

        q = session.query(
                cls,
                Experiment,
                Region,
                Source
            )\
            .join()\
            .filter(
                Experiment.uid == cls.experimentID,
                Region.uid == cls.regionID,
                Source.uid == cls.sourceID
            ).filter(
                cls.gene != None
            )

        if as_genomic_feature:

            feats = []

            # For each feature...
            for feat in q.all():
                feats.append(
                    cls.__as_genomic_feature(feat)
                )

            return feats

        return q.all()


#    @classmethod
#    def select_by_gene_tss(cls, session, gene,
#        tss, as_genomic_feature=False):
#        """
#        Query objects by gene TSS.
#        """
#
#        q = session.query(
#                cls,
#                Experiment,
#                Region,
#                Source
#            )\
#            .join()\
#            .filter(
#                Experiment.uid == cls.experimentID,
#                Region.uid == cls.regionID,
#                Source.uid == cls.sourceID
#            ).filter(
#                cls.gene == gene,
#                cls.tss == tss
#            )
#
#        if as_genomic_feature:
#
#            feats = []
#
#            # For each feature...
#            for feat in q.all():
#                feats.append(
#                    cls.__as_genomic_feature(feat)
#                )
#
#            return cls.__as_genomic_feature(
#                q.first()
#            )
#
#        return q.first()
#
#    @classmethod
#    def select_by_multiple_tss(cls, session, tss=[]):
#        """
#        Query objects by multiple TSSs. If no TSSs
#        are provided, return all objects. TSSs are
#        to be provided as a two-dimensional list in
#        the form: [[geneA, tss1], [geneA, tss2], ...]
#        """
#
#        # Initialize
#        ands = []
#
#        # For each gene, TSS pair...
#        for i, j in tss:
#            ands.append(
#                and_(
#                    cls.gene == i,
#                    cls.tss == j
#                )
#            )
#
#        q = session.query(cls).filter(or_(*ands))
#
#        return q.all()

    @classmethod
    def __as_genomic_feature(self, feat):

        # Define qualifiers
        qualifiers = {
            "uid": feat.TSS.uid,
            "regionID": feat.TSS.regionID,
            "gene": feat.TSS.gene,
            "tss": feat.TSS.tss,
            "sampleIDs": feat.TSS.sampleIDs,
            "avg_expression_levels": feat.TSS.avg_expression_levels,
            "experimentID": feat.TSS.experimentID,
            "sourceID": feat.TSS.sourceID,
            "experiment": feat.Experiment.name,
            "source" : feat.Source.name,            
        }

        if feat.TSS.gene:
            feat_id = "p%s@%s" % (
                feat.TSS.tss,
                feat.TSS.gene
            )
        else:
            feat_id = "p@%s:%s..%s,%s" % (
                feat.Region.chrom,
                int(feat.Region.start),
                int(feat.Region.end),
                feat.Region.strand
            )

        return GenomicFeature(
            feat.Region.chrom,
            int(feat.Region.start),
            int(feat.Region.end),
            strand = feat.Region.strand,
            feat_type = "TSS",
            feat_id = feat_id,
            qualifiers = qualifiers
        )

    def __repr__(self):

        return "<TSS(%s, %s, %s, %s, %s, %s, %s, %s)>" % \
            (
                "uid={}".format(self.uid),
                "regionID={}".format(self.regionID),
                "gene={}".format(self.gene),
                "tss={}".format(self.tss),
                "sampleIDs={}".format(self.sampleIDs),
                "avg_expression_levels={}".format(
                    self.avg_expression_levels
                ),
                "experimentID={}".format(self.experimentID),
                "sourceID={}".format(self.sourceID)
            )