from binning import (
    assign_bin,
    containing_bins,
    contained_bins
)
from sqlalchemy import (
    CheckConstraint,
    Column,
    Index,
    PrimaryKeyConstraint,
    String,
    ForeignKey,
    UniqueConstraint,
)
from sqlalchemy.dialects import mysql

from .base import Base
from .chrom import Chrom

class Region(Base):

    __tablename__ = "regions"

    uid = Column(
        "uid",
        mysql.INTEGER(unsigned=True)
    )

    bin = Column(
        "bin",
        mysql.SMALLINT(unsigned=True),
        nullable=False
    )

    chrom = Column(
        "chrom",
        String(5),
        ForeignKey("chroms.chrom"),
        nullable=False
    )

    start = Column(
        "start",
        mysql.INTEGER(unsigned=True),
        nullable=False
    )

    end = Column(
        "end",
        mysql.INTEGER(unsigned=True),
        nullable=False
    )

    strand = Column(
        "strand",
        mysql.CHAR(1)
    )

    __table_args__ = (
        PrimaryKeyConstraint(uid),
        UniqueConstraint(
            chrom,
            start,
            end,
            strand
        ),
        CheckConstraint("end > start"),
        Index("ix_bin_chrom", bin, chrom),
        {
            "mysql_engine": "MyISAM",
            "mysql_charset": "utf8"
        }
    )

    @classmethod
    def is_unique(cls, session, chrom, start, end,
        strand=None):

        q = session.query(cls).\
            filter(
                cls.chrom == chrom,
                cls.start == int(start),
                cls.end == int(end),
                cls.strand == strand
            )

        return len(q.all()) == 0

    @classmethod
    def select_unique(cls, session, chrom, start, end,
        strand=None):

        q = session.query(cls).\
            filter(
                cls.chrom == chrom,
                cls.start == int(start),
                cls.end == int(end),
                cls.strand == strand
            )

        return q.first()

    @classmethod
    def select_by_bin_range(cls, session, chrom,
        start, end, bins=[], compute_bins=False):
        """
        Query objects using the bin system to speed
        up range searches. If no bins are provided
        and compute_bins is set to True, then compute
        them. Otherwise, perform the query without
        using the bin system (EXTREMELY slow!).
        """

        if not bins and compute_bins:
            bins = cls._compute_bins(start, end)

        q = session.query(cls).\
            filter(
                cls.chrom == chrom,
                cls.end > start,
                cls.start < end
            )

        if bins:
            q = q.filter(cls.bin.in_(bins))

        return q.all()

    @classmethod
    def _compute_bins(cls, start, end):

        return list(set(
                containing_bins(start, end) +\
                contained_bins(start, end)
            )
        )

#    @classmethod
#    def select_by_exact_location(cls, session, chrom,
#        start, end):
#
#        bin = assign_bin(start, end)
#
#        q = session.query(cls).\
#            filter(
#                cls.bin == bin,
#                cls.chrom == chrom,
#                cls.start == start,
#                cls.end == end
#            )
#
#        return q.first()

    def __str__(self):

        return "{}\t{}\t{}\t{}\t{}".\
            format(
                self.bin,
                self.chrom,
                self.start,
                self.end,
                self.strand
            )

    def __repr__(self):

        return "<Region(%s, %s, %s, %s, %s, %s)>" % \
            (
                "uid={}".format(self.uid),
                "bin={}".format(self.bin),
                "chrom={}".format(self.chrom),
                "start={}".format(self.start),
                "end={}".format(self.end),
                "strand={}".format(self.strand)
            )