from binning import (
    containing_bins,
    contained_bins
)
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