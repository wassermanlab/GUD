from sqlalchemy import (
    Column, Index, PrimaryKeyConstraint, String, ForeignKey,
    UniqueConstraint, Integer, Float
)
from sqlalchemy.dialects import mysql
from GUD2.ORM.region import Region
from GUD2.ORM.source import Source
from GUD2.ORM.sample import Sample
from GUD2.ORM.experiment import Experiment
from GUD2.ORM.gene import Gene 
from GUD2.ORM.base import Base
from binning import containing_bins, contained_bins

class TSS(Base):

    __tablename__ = "tss"

    uid = Column("uid", mysql.INTEGER(unsigned=True))
    regionID = Column("regionID", Integer, ForeignKey("regions.uid"), nullable=False)
    sourceID = Column("sourceID", Integer, ForeignKey("sources.uid"), nullable=False)
    sampleID = Column("sampleID", Integer, ForeignKey("samples.uid"), nullable=False)
    strand = Column("strand", mysql.CHAR(1), nullable=False)
    experimentID = Column("experimentID", Integer, ForeignKey("experiments.uid"), nullable=False)
    gene = Column("gene", String(75), ForeignKey("genes.name2"))
    tss = Column("tss", mysql.INTEGER(unsigned=True))
    avg_tpm = Column("avg_tpm", Float, nullable=False)

    __table_args__ = (
        PrimaryKeyConstraint(uid),
        UniqueConstraint(regionID, sourceID, sampleID, experimentID),

        Index("ix_tss", regionID), ## query by bin range 
        Index("ix_tss_gene", gene, tss),

        {
            "mysql_engine": "MyISAM",
            "mysql_charset": "utf8"
        }
    )

    @classmethod
    def is_unique(cls, session, regionID, sourceID, sampleID, experimentID):

        q = session.query(cls).filter(
            cls.regionID == regionID,
            cls.sourceID == sourceID,
            cls.sampleID == sampleID,
            cls.experimentID == experimentID
        )

        return len(q.all()) == 0