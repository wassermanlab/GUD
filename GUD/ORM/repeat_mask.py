from binning import containing_bins, contained_bins 

from sqlalchemy import (
    Column, Date, Index, PrimaryKeyConstraint, String
)

from sqlalchemy.dialects import mysql
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.ext.hybrid import hybrid_property

Base = declarative_base()

class RepeatMask(Base):

    __tablename__ = "rmsk"

    bin = Column("bin", mysql.SMALLINT(unsigned=True), nullable=False)
    swScore = Column("swScore", mysql.INTEGER(unsigned=True), nullable=False)
#    milliDiv = Column("milliDiv", mysql.INTEGER(unsigned=True), nullable=False)
#    milliDel = Column("milliDel", mysql.INTEGER(unsigned=True), nullable=False)
#    milliIns = Column("milliIns", mysql.INTEGER(unsigned=True), nullable=False)
    genoName = Column("genoName", String(5), nullable=False)
    genoStart = Column("genoStart", mysql.INTEGER(unsigned=True), nullable=False)
    genoEnd = Column("genoEnd", mysql.INTEGER(unsigned=True), nullable=False)
#    genoLeft = Column("genoLeft", mysql.INTEGER(), nullable=False, default=0)
    strand = Column("strand", mysql.CHAR(1), nullable=False)
    repName = Column("repName", String(75), nullable=False)
    repClass = Column("repClass", String(75), nullable=False)
    repFamily = Column("repFamily", String(75), nullable=False)
#    repStart = Column("repStart", mysql.INTEGER(), nullable=False, default=0)
#    repEnd = Column("repEnd", mysql.INTEGER(), nullable=False, default=0)
#    repLeft = Column("repLeft", mysql.INTEGER(), nullable=False, default=0)
#    id = Column("id", mysql.CHAR(1), nullable=False)
    source_name = Column("source_name", String(25), nullable=False)
    date = Column("date", Date(), nullable=True)

    __table_args__ = (

        PrimaryKeyConstraint(
            genoName, genoStart, genoEnd, strand, repName,
            repClass, repFamily, source_name
        ),

        Index("ix_rmsk", bin, genoName),

        {
            "mysql_engine": "MyISAM",
            "mysql_charset": "utf8"
        }
    )

    # Create these properties to map columns to
    # standard attributes
    @hybrid_property
    def chrom(self):
        return self.genoName

    @chrom.setter
    def chrom(self, val):
        self.genoName = val

    @hybrid_property
    def start(self):
        return self.genoStart

    @start.setter
    def start(self, val):
        self.genoStart = val

    @hybrid_property
    def end(self):
        return self.genoEnd

    @end.setter
    def end(self, val):
        self.genoEnd = val

    @classmethod
    def select_by_bin_range(cls, session, chrom, start, end,
        repeat_classes=[], bins=[], compute_bins=False):
        """
        Query objects by chromosomal range using the binning system to
        speed up range searches. If bins are provided, use the given bins.
        If bins are NOT provided AND compute_bins is set to True, then
        compute the bins. Otherwise, perform the range query without the use
        of bins.
        """

        if not bins and compute_bins:
            bins = containing_bins(start, end) + contained_bins(start, end)

        q = session.query(cls).filter(
                cls.chrom == chrom, cls.end > start, cls.start < end)

        if bins:
            q = q.filter(cls.bin.in_(bins))

        if repeat_classes:
            q = q.filter(cls.repClass.in_(repeat_classes))

        return q.all()

    def __str__(self):
        return "{}\t{}\t{}\t{}|{}|{}\t{}\t{}".format(self.chrom,
            self.start, self.end, self.repName, self.repClass,
            self.repFamily, self.swScore, self.strand)

    def __repr__(self):
        return "<RepeatRegion(chrom={}, start={}, end={}, repName={}, strand={}, repClass={}, repFamily={}, score={})>".format(
            self.chrom, self.start, self.end, self.repName, self.strand,
            self.repClass, self.repFamily, self.swScore)