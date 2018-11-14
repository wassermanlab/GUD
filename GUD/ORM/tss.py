from binning import containing_bins, contained_bins 
from Bio.SeqFeature import FeatureLocation

from sqlalchemy import (
    Column, Date, Enum, Index, Integer,
    PrimaryKeyConstraint, String
)

from sqlalchemy.dialects import mysql
from sqlalchemy.ext.declarative import declarative_base

Base = declarative_base()

class TSS(Base):

    __tablename__ = "tss"

    bin = Column("bin", mysql.SMALLINT(unsigned=True), nullable=False)
    gene = Column("gene", String(75))
    id = Column("id", Integer)
    chrom = Column("chrom", String(5), nullable=False)
    start = Column("start", mysql.INTEGER(unsigned=True), nullable=False)
    end = Column("end", mysql.INTEGER(unsigned=True), nullable=False)
    strand = Column("strand", mysql.CHAR(1), nullable=False)
    tpm = Column("tpm", mysql.LONGBLOB, nullable=False)
    source_name = Column("source_name", String(25), nullable=False)
    date = Column("date", Date(), nullable=True)

    __table_args__ = (

        PrimaryKeyConstraint(
            gene, id, chrom, start, end, strand, source_name
        ),

        Index("ix_tss", bin, chrom),
        Index("ix_tss_name", name),

        {
            "mysql_engine": "MyISAM",
            "mysql_charset": "utf8"
        }
    )

#    def __str__(self):
#        return "{}\t{}\t{}\tp{}{}".format(self.chrom, self.start, self.end,
#            self.relative_expression)
#
#    def __repr__(self):
#        return "<Expression(cage_peak_id={}, sample_id={}, tag_count={}, tpm={}, relative_expression={})>".format(self.cage_peak_id, self.sample_id, self.tag_count, self.tpm, self.relative_expression)