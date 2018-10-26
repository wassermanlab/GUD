from ..utils.bin_range import BinRange

from sqlalchemy import (
    Column, Date, Index, PrimaryKeyConstraint, String
)

from sqlalchemy.dialects import mysql
from sqlalchemy.ext.declarative import declarative_base

Base = declarative_base()

class Tad(Base):

    __tablename__ = "tad"

    bin = Column("bin", mysql.SMALLINT(unsigned=True), nullable=False)
    chrom = Column("chrom", String(5), nullable=False)
    start = Column("start", mysql.INTEGER(unsigned=True), nullable=False)
    end = Column("end", mysql.INTEGER(unsigned=True), nullable=False)
    hierarchical_level = Column("hierarchical_level", mysql.SMALLINT(unsigned=True),
        nullable=False)
    cell_or_tissue = Column("cell_or_tissue", String(225), nullable=False)
    experiment_type = Column("experiment_type", String(25), nullable=False)
    restriction_enzyme = Column("restriction_enzyme", String(25), nullable=False)
    source_name = Column("source_name", String(25), nullable=False)
    date = Column("date", Date(), nullable=True)

    __table_args__ = (

        PrimaryKeyConstraint(
            bin, chrom, start, end, hierarchical_level, cell_or_tissue,
            experiment_type, restriction_enzyme, source_name
        ),

        Index("ix_tad", bin, chrom),

        {
            "mysql_engine": "MyISAM",
            "mysql_charset": "utf8"
        }
    )

    @classmethod
    def select_by_bin_range(cls, session, chrom, start, end,
        sample=[], bins=[], compute_bins=False):
        """
        Query objects by chromosomal range using the binning system to
        speed up range searches. If bins are provided, use the given bins.
        If bins are NOT provided AND compute_bins is set to True, then
        compute the bins. Otherwise, perform the range query without the use
        of bins.
        """

        if not bins and compute_bins:
            bins = BinRange().allBinsInRange(start, end)

        q = session.query(cls).filter(
                cls.chrom == chrom, cls.end > start, cls.start < end)

        if bins:
            q = q.filter(cls.bin.in_(bins))

        if sample:
            q = q.filter(cls.cell_or_tissue.in_(sample))

        return q.all()

    @classmethod
    def select_encompasing_range(cls, session, chrom, start, end,
        sample=[], bins=[], compute_bins=False):
        """
        Query objects which encompass the given range, e.g. an object
        which fully contains the given feature.
        """

        if not bins and compute_bins:
            bins = BinRange().allBinsInRange(start, end)

        q = session.query(cls).filter(
                cls.chrom == chrom, cls.start <= start, cls.end >= end)

        if bins:
            q = q.filter(cls.bin.in_(bins))

        if sample:
            q = q.filter(cls.cell_or_tissue.in_(sample))

        return q.all()

    @classmethod
    def select_overlapping_range(cls, session, chrom, start, end,
        sample=[], bins=[], compute_bins=False):
        """
        Query objects which overlap the given range, e.g. an object
        which partially overlaps the given feature.
        """

        if not bins and compute_bins:
            bins = BinRange().allBinsInRange(start, end)

        q = session.query(cls).filter(
                cls.chrom == chrom, cls.start < end, cls.end > start)

        if bins:
            q = q.filter(cls.bin.in_(bins))

        if sample:
            q = q.filter(cls.cell_or_tissue.in_(sample))

        return q.all()

    @classmethod
    def feature_exists(cls, session, chrom, start, end,
        hierarchical_level, cell_or_tissue, experiment_type,
        restriction_enzyme, source_name): 
        """
        Returns whether a feature exists in the database.
        """

        q = session.query(cls).filter(
                cls.chrom == chrom,
                cls.start == start,
                cls.end == end,
                cls.hierarchical_level == hierarchical_level,
                cls.cell_or_tissue == cell_or_tissue,
                cls.experiment_type == experiment_type,
                cls.restriction_enzyme == restriction_enzyme,
                cls.source_name == source_name
            )

        return session.query(q.exists()).scalar()

    def __str__(self):
        return "{}\t{}\t{}\t{}|{}|{}|{}\t{}".format(self.chrom, self.start,
            self.end, self.experiment_type, self.cell_or_tissue,
            self.source_name, self.restriction_enzyme, self.hierarchical_level)

    def __repr__(self):
        return "<Tad(chrom={}, start={}, end={}, level={}, sample={}, experiment={}, enzyme={}, source={}, date={})>".format(
            self.chrom, self.start, self.end, self.hierarchical_level,
            self.cell_or_tissue, self.experiment_type, self.restriction_enzyme,
            self.source_name, self.date)
