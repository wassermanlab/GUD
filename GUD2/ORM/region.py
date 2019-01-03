from binning import assign_bin, containing_bins, contained_bins

from sqlalchemy import (
    Column, Index, PrimaryKeyConstraint, String,
    ForeignKey, UniqueConstraint, CheckConstraint
)

from sqlalchemy.dialects import mysql

from GUD2.ORM.base import Base

class Region(Base):

    __tablename__ = "regions"

    uid = Column("uid", mysql.INTEGER(unsigned=True))

    bin = Column("bin", mysql.SMALLINT(unsigned=True),
        nullable=False)

    chrom = Column("chrom", String(5),
        ForeignKey("chroms.chrom"), nullable=False)

    start = Column("start",
        mysql.INTEGER(unsigned=True), nullable=False)

    end = Column("end", mysql.INTEGER(unsigned=True),
        nullable=False)

    __table_args__ = (

        PrimaryKeyConstraint(uid),
        UniqueConstraint(chrom, start, end),
        CheckConstraint("end > start"),
        Index("ix_region", bin, chrom),

        {
            "mysql_engine": "MyISAM",
            "mysql_charset": "utf8"
        }
    )

    @classmethod
    def select_by_bin_range(cls, session, chrom, start,
        end, bins=[], compute_bins=False, as_list=True):
        """
        Query objects by chromosomal range using the
        binning system to speed up searches. If \"bins\"
        are provided, use them. If \"bins\" are NOT
        provided AND \"compute_bins\" is \"True\",
        compute them. Otherwise, query objects without
        using the binning system.
        """

        if not bins and compute_bins:
            bins = containing_bins(start, end)
            bins += contained_bins(start, end)

        q = session.query(cls).filter(
            cls.chrom == chrom,
            cls.end > start,
            cls.start < end
        )

        if bins:
            q = q.filter(cls.bin.in_((list(bins))))

        if return_list:
            return q.all()

        return q

    @classmethod
    def select_unique(cls, session, chrom, start, end):

        q = session.query(cls).filter(
            cls.bin == assign_bin(start, end),
            cls.chrom == chrom,
            cls.start == start,
            cls.end == end)

        return q.first()

    def __str__(self):
        return "{}\t{}\t{}\t{}".format(
            self.bin,
            self.chrom,
            self.start,
            self.end,
        )

    def __repr__(self):
        return "<Region(uid={}, bin={}, chrom={}, start={}, end={})>".format(
            self.uid,
            self.bin,
            self.chrom,
            self.start,
            self.end
        )