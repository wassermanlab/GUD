from sqlalchemy import (
    Column, Index, PrimaryKeyConstraint, String, ForeignKey,
    UniqueConstraint, CheckConstraint, Integer, Float
)
from sqlalchemy.dialects import mysql
from GUD2.ORM.region import Region
from GUD2.ORM.source import Source
from GUD2.ORM.sample import Sample
from GUD2.ORM.experiment import Experiment
from GUD2.ORM.gene import Gene ## todo
from GUD2.ORM.base import Base
from binning import containing_bins, contained_bins

class TSS(Base):

    __tablename__ = "tss"

    uid = Column("uid", mysql.INTEGER(unsigned=True))
    regionID = Column("regionID", Integer, ForeignKey('regions.uid'), nullable=False)
    sourceID = Column("sourceID", Integer, ForeignKey('sources.uid'), nullable=False)
    sampleID = Column("sampleID", Integer, ForeignKey('samples.uid'), nullable=False)
    experimentID = Column("experimentID", Integer, ForeignKey('experiments.uid'), nullable=False)
    gene = Column("gene", String(75), ForeignKey('genes.name2'))
    tss = Column("tss", mysql.INTEGER(unsigned=True))
    avg_tpm = Column("avg_tpm", Float, nullable=False)
    strand = Column("strand", mysql.CHAR(1), nullable=False)

    __table_args__ = (
        PrimaryKeyConstraint(uid),
        UniqueConstraint(uid, regionID, sourceID, sampleID, experimentID),

        Index("ix_tss", regionID), ## query by bin range 
        Index("ix_tss_gene", gene, tss),

        {
            "mysql_engine": "MyISAM",
            "mysql_charset": "utf8"
        }
    )

    def __str__(self):
        return "{}\t{}".format(self.motif, self.pathogenicity)

    def __repr__(self):
        return "<ShortTandemRepeat(uid={}, regionID={}, sourceID={}, motif={}, pathogencity={})>".format(
            self.uid, self.regionID, self.sourceID, self.motif, self.pathogenicity)
