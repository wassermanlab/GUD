from binning import containing_bins, contained_bins 
from Bio.SeqFeature import FeatureLocation

from sqlalchemy import (
    Column, Date, Enum, Float, Index,
    Integer, PrimaryKeyConstraint, String, types
)

from sqlalchemy.dialects import mysql
from sqlalchemy.ext.declarative import declarative_base

Base = declarative_base()

class TSS(Base):

    __tablename__ = "tss"

    bin = Column("bin", mysql.SMALLINT(unsigned=True), nullable=False)
    gene = Column("gene", String(75))
    tss = Column("tss", Integer)
    chrom = Column("chrom", String(5), nullable=False)
    start = Column("start", mysql.INTEGER(unsigned=True), nullable=False)
    end = Column("end", mysql.INTEGER(unsigned=True), nullable=False)
    strand = Column("strand", mysql.CHAR(1), nullable=False)
    cell_or_tissue = Column("cell_or_tissue", String(225), nullable=False)
    replicate = Column("replicate", Integer, nullable=False)
    tpm = Column("tpm", Float, nullable=False)
    experiment_type = Column("experiment_type", String(25), nullable=False)
    source_name = Column("source_name", String(25), nullable=False)
    date = Column("date", Date(), nullable=True)

    __table_args__ = (

        PrimaryKeyConstraint(
            chrom, start, end, strand, cell_or_tissue, replicate,
            experiment_type, source_name
        ),

        Index("ix_tss", bin, chrom),
        Index("ix_tss_gene", gene),
        Index("ix_tss_gene_tss", gene, tss),
        Index("ix_tss_cell_or_tissue", cell_or_tissue),

        {
            "mysql_engine": "MyISAM",
            "mysql_charset": "utf8"
        }
    )

    def __str__(self):
        return "{}\t{}\t{}\t{} ({})\t{}\t{}\t{} ({})".format(self.chrom,
            self.start, self.end, self.gene, self.tss, self.strand, self.tpm,
            self.cell_or_tissue, self.replicate)

    def __repr__(self):
        return "<TSS(gene={}, tss={}, chrom={}, start={}, end={}, strand={}, sample={}, replicate={}, tpm={}, experiment={}, source={}, date={})>".format(
            self.gene, self.tss, self.chrom, self.start, self.end, self.strand, self.cell_or_tissue, self.replicate, self.tpm, self.experiment, self.source_name, self.date)