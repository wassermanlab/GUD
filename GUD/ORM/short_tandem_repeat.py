from binning import containing_bins, contained_bins 
from Bio.SeqFeature import FeatureLocation

from sqlalchemy import (
    Column, Index, Integer,
    PrimaryKeyConstraint, String
)

from datetime import date
from sqlalchemy.dialects import mysql
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.ext.hybrid import hybrid_property

Base = declarative_base()


class ShortTandemRepeat(Base):
    __tablename__ = 'short_tandem_repeat'

    bin = Column("bin", mysql.SMALLINT(unsigned=True), nullable=False)
    chrom = Column("chrom", String(30), nullable=False)
    start = Column("start", mysql.INTEGER(unsigned=True), nullable=False)
    end = Column("end", mysql.INTEGER(unsigned=True), nullable=False)
    length = Column("length", mysql.INTEGER(unsigned=True), nullable=False)
    motif = Column("motif", String(5), nullable=False)
    pathogenicity = Column("pathogenicity", mysql.INTEGER(
        unsigned=False), nullable=False)
    date = Column("date", Date(), nullable=True)

    __table_args__ = (

        PrimaryKeyConstraint(
            chrom, start, end
        ),

        Index("ix_str", bin, chrom),
        Index("ix_str_pathogenicity", pathogenicity),

        {
            "mysql_engine": "MyISAM",
            "mysql_charset": "utf8"
        }
    )

    def __str__(self):
        return "{}\t{}\t{}\t{}\t{}".format(
            self.chrom, self.start, self.end, self.length, self.motif)

    def __repr__(self):
        return "<ShortTandemRepeat(bin={}, chrom={}, start={}, end={}, length={}, motif={}, pathogenicity={})>".format(
            self.bin, self.chrom, self.start, self.end, self.length, self.motif, self.pathogenicity)
