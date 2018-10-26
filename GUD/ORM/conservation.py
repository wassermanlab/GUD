from binning import containing_bins

from sqlalchemy import (
    Column, Date, Index, PrimaryKeyConstraint, String
)

from sqlalchemy.dialects import mysql
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.ext.hybrid import hybrid_property

Base = declarative_base()

class Conservation(Base):
        
    __tablename__ = "conservation"

    bin = Column("bin", mysql.SMALLINT(unsigned=True), nullable=False)
    chrom = Column("chrom", String(5), nullable=False)
    chromStart = Column("chromStart", mysql.INTEGER(unsigned=True), nullable=False)
    chromEnd = Column("chromEnd", mysql.INTEGER(unsigned=True), nullable=False)
#    extFile = Column("extFile", mysql.INTEGER(), nullable=False, default=0)
#    offset = Column("offset", mysql.INTEGER(), nullable=False, default=0)
    score = Column("score", mysql.FLOAT(unsigned=True), nullable=False)
    source_name = Column("source_name", String(25), nullable=False)
    date = Column("date", Date(), nullable=True) 

    __table_args__ = (

        PrimaryKeyConstraint(
            chrom, chromStart, chromEnd, source_name
        ),

        Index("ix_conservation", bin, chrom),
        
        {
            "mysql_engine": "MyISAM",
            "mysql_charset": "utf8"
        }
    )

    # Create these properties to map columns to
    # standard attributes
    @hybrid_property
    def start(self):
        return self.chromStart

    @start.setter
    def start(self, val):
        self.chromStart = val

    @hybrid_property
    def end(self):
        return self.chromEnd

    @end.setter
    def end(self, val):
        self.chromEnd = val

    @classmethod
    def select_by_bin_range(cls, session, chrom, start, end, bins=[],
        compute_bins=False): 
        """
        Query objects by chromosomal range using the binning system to
        speed up range searches. If bins are provided, use the given bins.
        If bins are NOT provided AND compute_bins is set to True, then
        compute the bins. Otherwise, perform the range query without the use
        of bins.
        """

        if not bins and compute_bins:
            bins = containing_bins(start, end)

        q = session.query(cls).filter(
                cls.chrom == chrom, cls.end > start, cls.start < end)

        if bins:
            q = q.filter(cls.bin.in_(bins))

        return q.all()

    def __str__(self):
        return "{}\t{}\t{}\t.\t{}".format(self.chrom, self.start,
            self.end, self.score)

    def __repr__(self):
        return "<Conservation(chrom={}, start={}, end={}, score={}, source={}, date={})>".format(
            self.chrom, self.start, self.end, self.score, self.source_name,
            self.date)
