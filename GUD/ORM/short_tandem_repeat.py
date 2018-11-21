from binning import containing_bins, contained_bins 
from Bio.SeqFeature import FeatureLocation

from sqlalchemy import (
    Column, Date, Index, Integer,
    PrimaryKeyConstraint, String
)

from datetime import date
from sqlalchemy.dialects import mysql
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.ext.hybrid import hybrid_property

Base = declarative_base()

class ShortTandemRepeat(Base):

    __tablename__ = "short_tandem_repeat"

    bin = Column("bin", mysql.SMALLINT(unsigned=True), nullable=False)
    chrom = Column("chrom", String(5), nullable=False)
    start = Column("start", mysql.INTEGER(unsigned=True), nullable=False)
    end = Column("end", mysql.INTEGER(unsigned=True), nullable=False)
    strand = Column("strand", mysql.CHAR(1), nullable=False)
    motif = Column("motif", String(5), nullable=False)
    pathogenicity = Column("pathogenicity", mysql.INTEGER(
        unsigned=False), nullable=False)
    source_name = Column("source_name", String(25), nullable=False)
    date = Column("date", Date(), nullable=True)

    __table_args__ = (

        PrimaryKeyConstraint(
            chrom, start, end, strand, source_name
        ),

        Index("ix_str", bin, chrom),
        Index("ix_str_pathogenicity", pathogenicity),

        {
            "mysql_engine": "MyISAM",
            "mysql_charset": "utf8"
        }
    )

    @classmethod
    def select_by_bin_range(cls, session, chrom, start, end,
        bins=[], compute_bins=False):
        """
        Query objects by chromosomal range using the binning system to
        speed up range searches. If bins are provided, use the given bins.
        If bins are NOT provided AND compute_bins is set to True, then
        compute the bins. Otherwise, perform the range query without the use
        of bins.
        """

        if not bins and compute_bins:
            bins = set(containing_bins(start, end) + contained_bins(start, end))

        q = session.query(cls).filter(
                cls.chrom == chrom, cls.start >= start, cls.start < end)

        if bins:
            q = q.filter(cls.bin.in_((list(bins))))

        return q.all()
    
    @classmethod
    def select_by_motif(cls, motif, compute_rotations=False):
        """
        Query objects by motif without, calculate all the rotations of a motif
        if asked
        """
        if compute_rotations is True: #comput the rotations
            motif_len = len(motif)
            motif_temp = motif 
            motifs = []
            for i in range(motif_len):
                motif_temp = motif_temp[1:motif_len] +motif_temp[0]
                motifs.append(motif_temp)    
        else:    
            motifs = [motif]

        q = session.query(cls).filter(cls.motif.in_(motifs))

        return q.all()    

    def __str__(self):
        return "{}\t{}\t{}\t{}\t{}".format(
            self.chrom, self.start, self.end, self.length, self.motif)

    def __repr__(self):
        return "<ShortTandemRepeat(bin={}, chrom={}, start={}, end={}, length={}, motif={}, pathogenicity={})>".format(
            self.bin, self.chrom, self.start, self.end, self.length, self.motif, self.pathogenicity)