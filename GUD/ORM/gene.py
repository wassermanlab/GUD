from binning import containing_bins, contained_bins 
from Bio.SeqFeature import FeatureLocation

from sqlalchemy import (
    Column, Date, Enum, Index, Integer,
    PrimaryKeyConstraint, String
)

from sqlalchemy.dialects import mysql
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.ext.hybrid import hybrid_property

Base = declarative_base()

class Gene(Base):
    
    __tablename__ = "gene"

    bin = Column("bin", mysql.SMALLINT(unsigned=True), nullable=False)
    name = Column("name", String(75), nullable=False)
    chrom = Column("chrom", String(5), nullable=False)
    strand = Column("strand", mysql.CHAR(1), nullable=False)
    txStart = Column("txStart", mysql.INTEGER(unsigned=True), nullable=False)
    txEnd = Column("txEnd", mysql.INTEGER(unsigned=True), nullable=False)
    cdsStart = Column("cdsStart", mysql.INTEGER(unsigned=True), nullable=False)
    cdsEnd = Column("cdsEnd", mysql.INTEGER(unsigned=True), nullable=False)
#    exonCounts = Column("exonCounts", mysql.INTEGER(unsigned=True), nullable=False)
    exonStarts = Column("exonStarts", mysql.LONGBLOB, nullable=False)
    exonEnds = Column("exonEnds", mysql.LONGBLOB, nullable=False)
#    score = Column("score", Integer)
    name2 = Column("name2", String(75), nullable=False)
#    cdsStartStat = Column("cdsStartStat",
#        Enum("none", "unk", "incmpl", "cmpl"), nullable=False)
#    cdsEndStat = Column("cdsEndStat",
#        Enum("none", "unk", "incmpl", "cmpl"), nullable=False)
#    exonFrames = Column("exonFrames", mysql.LONGBLOB, nullable=False)
    source_name = Column("source_name", String(25), nullable=False)
    date = Column("date", Date(), nullable=True)

    _exons = []
    _coding_exons = []

    __table_args__ = (

        PrimaryKeyConstraint(
            name, chrom, strand, txStart, txEnd, source_name
        ),

        Index("ix_gene", bin, chrom),
        Index("ix_gene_name", name),
        Index("ix_gene_name2", name2),

        {
            "mysql_engine": "MyISAM",
            "mysql_charset": "utf8"
        }
    )

    # Create these properties to map columns to
    # standard attributes
    @hybrid_property
    def start(self):
        return self.txStart

    @start.setter
    def start(self, val):
        self.txStart = val

    @hybrid_property
    def end(self):
        return self.txEnd

    @end.setter
    def end(self, val):
        self.txEnd = val

    @property
    def exons(self):
        if not self._exons:
            self.compute_exons()

        return self._exons

    @property
    def coding_exons(self):
        if not self._coding_exons:
            self.compute_coding_exons()

        return self._coding_exons

    @classmethod
    def select_all_gene_symbols(cls, session):
        """
        Query all gene objects in the database and return their gene
        symbols (name2 field).
        """

        all_genes = session.query(cls).all()

        symbols = []
        for g in all_genes:
            if g.name2:
                symbols.append(g.name2)

        symbols.sort()

        return symbols

    @classmethod
    def select_by_bin_range(cls, session, chrom, start, end, bins=[],
        compute_bins=False):
        """
        Query objects by chromosomal range using the binning system to
        speed up range searches. If bins are provided use the given bins.
        If bins are NOT provided AND compute_bins is set to True, then
        compute the bins. Otherwise perform the range query without the use
        of bins.
        """

        if not bins and compute_bins:
            bins = containing_bins(start, end) + contained_bins(start, end)

        q = session.query(cls).filter(
                cls.chrom == chrom, cls.end > start, cls.start < end)

        if bins:
            q = q.filter(cls.bin.in_(bins))

        return q.all()

    @classmethod
    def select_overlapping_range(cls, session, chrom, start, end, tissue=[]):
        """
        Query genes which overlap the given range, e.g. a gene which
        at least partially overlaps the given feature.
        """

        q = session.query(cls).filter(
                cls.chrom == chrom, cls.start < end, cls.end > start)

        if tissue:
            q = q.filter(cls.tissue.in_(tissue))

        return q.all()

    @classmethod
    def select_by_name(cls, session, name):
        """
        Query refGene objects by common name. If no name is provided,
        query all genes.
        """

        q = session.query(cls)

        if name:
            q = q.filter(cls.name2 == name)

        return q.all()

    @classmethod
    def select_by_names(cls, session, names=[]):
        """
        Query refGene objects by list of common names. If no names are
        provided, query all genes.
        """

        q = session.query(cls)

        if names:
            q = q.filter(cls.name2.in_(names))

        return q.all()

    def __str__(self):
        return "{}\t{}\t{}\t{}\t0\t{}".format(self.chrom, self.txStart, self.txEnd, self.name2, self.strand)

    def __repr__(self):
        return "<Gene(name={}, name2={}, chrom={}, txStart={}, txEnd={}, strand={}, source={}, date={})>".format(
            self.name, self.name2, self.chrom, self.txStart, self.txEnd, self.strand, self.source_name, self.date)