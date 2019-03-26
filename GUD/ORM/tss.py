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
        UniqueConstraint(
            regionID,
            experimentID,
            sourceID
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
        experimentID):

        q = session.query(cls).\
            filter(
                cls.regionID == regionID,
                cls.sourceID == sourceID,
                cls.experimentID == experimentID
            )

        return len(q.all()) == 0


    @classmethod
    def select_unique(cls, session, regionID,
        sourceID, experimentID):

        q = session.query(cls).\
            filter(
                cls.regionID == regionID,
                cls.sourceID == sourceID,
                cls.experimentID == experimentID
            )

        return q.first()

    @classmethod
    def select_by_uid(cls, session, uid):

        q = session.query(cls).\
            filter(
                cls.uid == uid
            )

        return q.first()

    @classmethod
    def select_by_tss(cls, session, gene, tss):
        """
        Query objects by TSS (i.e. gene + tss).
        """
    
        q = session.query(cls).\
            filter(
                cls.gene == gene,
                cls.tss == tss
            )

        return q.first()

    @classmethod
    def select_by_multiple_tss(cls, session, tss=[]):
        """
        Query objects by multiple TSSs. If no TSSs
        are provided, return all objects. TSSs are
        to be provided as a two-dimensional list in
        the form: [[geneA, tss1], [geneA, tss2], ...]
        """

        # Initialize
        ands = []

        # For each gene, TSS pair...
        for i, j in tss:
            ands.append(
                and_(
                    cls.gene == i,
                    cls.tss == j
                )
            )

        q = session.query(cls).filter(or_(*ands))

        return q.all()

    def __as_genomic_feature(feat):

        # Initialize
        exonStarts = []
        exonEnds = []

        # For each exon start...
        for i in str(feat.Gene.exonStarts).split(","):
            if i.isdigit():
                exonStarts.append(int(i))
        
        # For each exon end...
        for i in str(feat.Gene.exonEnds).split(","):
            if i.isdigit():
                exonStarts.append(int(i))

        # Define qualifiers
        qualifiers = {
            "name": feat.Gene.name,
            "cdsStart": int(feat.Gene.cdsStart),
            "cdsEnd": int(feat.Gene.cdsEnd),
            "exonStarts": exonStarts,
            "exonEnds": exonEnds,
            "source" : feat.Source.name,            
        }

        return GenomicFeature(
            feat.Region.chrom,
            int(feat.Region.start),
            int(feat.Region.end),
            strand = feat.Region.strand,
            feat_type = "Gene",
            feat_id = feat.Gene.name2,
            qualifiers = qualifiers
        )

#    def __repr__(self):
#
#        return "<TSS(%s, %s, %s, %s)>" % \
#            (
#                "uid={}".format(self.uid),
#                "regionID={}".format(self.regionID),
#                "name={}".format(self.name),
#                "name2={}".format(self.name2),
#                "sourceID={}".format(self.sourceID),
#            )