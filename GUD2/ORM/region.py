from sqlalchemy import (
    Column, Index, PrimaryKeyConstraint, String, ForeignKey,
    UniqueConstraint, CheckConstraint
)
from sqlalchemy.dialects import mysql
from binning import assign_bin, containing_bins, contained_bins
from GUD2.ORM.chrom import Chrom
from GUD2.ORM.base import Base

class Region(Base):

    __tablename__ = "regions"

    uid = Column("uid", mysql.INTEGER(unsigned=True))
    bin = Column("bin", mysql.SMALLINT(unsigned=True), nullable=False)
    chrom = Column("chrom", String(5), ForeignKey('chroms.chrom'), nullable=False)
    start = Column("start", mysql.INTEGER(unsigned=True), nullable=False)
    end = Column("end", mysql.INTEGER(unsigned=True), nullable=False)

    __table_args__ = (
        PrimaryKeyConstraint(uid),
        UniqueConstraint(chrom, start, end),
        CheckConstraint('end > start'),

        Index("ix_region", bin, chrom),

        {
            "mysql_engine": "MyISAM",
            "mysql_charset": "utf8"
        }
    )

    @classmethod
    def select_by_bin_range(cls, session, chrom, start, end,
                            bins=[], compute_bins=False, return_list=True):
        """
        Query objects by chromosomal range using the binning system to
        speed up range searches. If bins are provided, use the given bins.
        If bins are NOT provided AND compute_bins is set to True, then
        compute the bins. Otherwise, perform the range query without the use
        of bins.
        """

        if not bins and compute_bins:
            bins = set(containing_bins(start, end) +
                       contained_bins(start, end))

        q = session.query(cls).filter(cls.chrom == chrom, cls.end > start,
                                      cls.start < end)

        if bins:
            q = q.filter(cls.bin.in_((list(bins))))

        if return_list:
            return q.all()

        return q

    @classmethod
    def select_by_exact_location(cls, session, chrom, start, end):

        bin = assign_bin(start, end)

        q = session.query(cls).filter(cls.bin == bin,
            cls.chrom == chrom, cls.start == start, cls.end == end)

        return q.first()

    def __str__(self):
        return "{}\t{}\t{}\t{}".format(self.bin, self.chrom,
            self.start, self.end)

    def __repr__(self):
        return "<Chrom(uid={}, bin={}, chrom={}, start={}, end={})>".format(
            self.uid, self.bin, self.chrom, self.start, self.end)