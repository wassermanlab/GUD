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

    uid = Column("uid", mysql.INTEGER(unsigned=True))
    regionID = Column("regionID", Integer, ForeignKey("regions.uid"), nullable=False)
    sourceID = Column("sourceID", Integer, ForeignKey("sources.uid"), nullable=False)
    experimentID = Column("experimentID", Integer, ForeignKey("experiments.uid"), nullable=False)
    gene = Column("gene", String(75), ForeignKey("genes.name2"))
    tss = Column("tss", mysql.INTEGER(unsigned=True))
    strand = Column("strand", mysql.CHAR(1), nullable=False)
    sampleIDs = Column("sampleIDs", mysql.LONGBLOB, nullable=False)
    avg_expression_levels = Column("avg_expression_levels", mysql.LONGBLOB, nullable=False)

    __table_args__ = (
        PrimaryKeyConstraint(uid),
        UniqueConstraint(
            regionID,
            sourceID,
            experimentID,
            gene,
            tss
        ),
        Index("ix_tss", regionID), # query by bin range 
        Index("ix_tss_gene", gene, tss),
        {
            "mysql_engine": "MyISAM",
            "mysql_charset": "utf8"
        }
    )

    @classmethod
    def is_unique(cls, session, regionID, sourceID,
        experimentID, gene, tss):

        q = session.query(cls).filter(
            cls.regionID == regionID,
            cls.sourceID == sourceID,
            cls.experimentID == experimentID,
            cls.gene == gene,
            cls.tss == tss
        )

        return len(q.all()) == 0

    @classmethod
    def select_by_exact_tss(cls, session, regionID,
        sourceID, experimentID, gene, tss):

        q = session.query(cls).filter(
            cls.regionID == regionID,
            cls.sourceID == sourceID,
            cls.experimentID == experimentID,
            cls.gene == gene,
            cls.tss == tss
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
        are provided, return all objects. Provide TSSs
        as a two-dimensinal list (i.e. [[gene, tss], ...]).
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