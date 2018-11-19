from binning import containing_bins, contained_bins 
from Bio.SeqFeature import FeatureLocation

from sqlalchemy import (
    Column, Date, Enum, Float, Index,
    Integer, PrimaryKeyConstraint, String, types
)

from sqlalchemy.dialects import mysql
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.sql import func

Base = declarative_base()

class TSS(Base):

    __tablename__ = "tss"

    bin = Column("bin", mysql.SMALLINT(unsigned=True), nullable=False)
    gene = Column("gene", String(75))
    tss = Column("tss", Integer)
    chrom = Column("chrom", String(5), nullable=False)
    start = Column("start", mysql.INTEGER(unsigned=True), nullable=False)
    end = Column("end", mysql.INTEGER(unsigned=True), nullable=False)
    strand = Column("strand", mysql.CHAR(1), nullable=False)
    cell_or_tissue = Column("cell_or_tissue", String(225), nullable=False)
    replicate = Column("replicate", Integer, nullable=False)
    tpm = Column("tpm", Float, nullable=False)
    percent_tpm = Column("percent_tpm", Float, nullable=False)
    experiment_type = Column("experiment_type", String(25), nullable=False)
    source_name = Column("source_name", String(25), nullable=False)
    date = Column("date", Date(), nullable=True)

    __table_args__ = (

        PrimaryKeyConstraint(
            chrom, start, end, strand, cell_or_tissue, replicate,
            experiment_type, source_name
        ),

        Index("ix_tss", bin, chrom),
        Index("ix_tss_gene", gene),
        Index("ix_tss_gene_tss", gene, tss),
        Index("ix_tss_cell_or_tissue", cell_or_tissue),

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
            bins = set(containing_bins(start, end) + contained_bins(start, end))

        q = session.query(cls).filter(
                cls.chrom == chrom, cls.end > start, cls.start < end)

        if bins:
            q = q.filter(cls.bin.in_((list(bins))))

        if sample:
            q = q.filter(cls.cell_or_tissue.in_(sample))

        return q.all()

    @classmethod
    def select_by_gene(cls, session, gene, sample=[]):
        """
        Query objects by gene.
        """

        q = session.query(cls).filter(cls.gene == gene)

        if sample:
            q = q.filter(cls.cell_or_tissue.in_(sample))

        return q.all()

    @classmethod
    def select_by_tss(cls, session, gene, tss, sample=[]):
        """
        Query objects by gene tss.
        """

        q = session.query(cls).filter(cls.gene == gene, cls.tss == tss)

        if sample:
            q = q.filter(cls.cell_or_tissue.in_(sample))

        return q.all()

    @classmethod
    def select_by_sample(cls, session, sample=[], min_tpm=0.0):
        """
        Query objects by sample with a minimum tpm.
        """

        q = session.query(func.avg(cls.tpm)).group_by(
            cls.chrom, cls.start, cls.end, cls.strand)

        print(q.all())
        exit(0)

        if sample:
            q = q.filter(cls.cell_or_tissue.in_(sample))

        q = q.query()

        return q.all()


        station_data = dbsession.query(func.hour(cls.last_update),
                                        func.avg(cls.available_bikes),
                                        func.avg(cls.available_bike_stands)) \
            .filter(cls.station_id == station_id,
                    func.weekday(cls.last_update) == weekday) \
            .group_by(func.hour(cls.last_update)) \
            .all()


    @classmethod
    def feature_exists(cls, session, chrom, start, end, strand,
        cell_or_tissue, replicate, experiment_type, source_name): 
        """
        Returns whether a feature exists in the database.
        """

        q = session.query(cls).filter(
                cls.chrom == chrom,
                cls.start == start,
                cls.end == end,
                cls.strand == strand,
                cls.cell_or_tissue == cell_or_tissue,
                cls.replicate == replicate,
                cls.experiment_type == experiment_type,
                cls.source_name == source_name
            )

        return session.query(q.exists()).scalar()

    def __str__(self):
        return "{}\t{}\t{}\t{} ({})\t{}\t{} ({}%)\t{} ({})".format(self.chrom,
            self.start, self.end, self.gene, self.tss, self.strand, self.tpm,
            self.percent_tpm, self.cell_or_tissue, self.replicate)

    def __repr__(self):
        return "<TSS(gene={}, tss={}, chrom={}, start={}, end={}, strand={}, sample={}, replicate={}, tpm={}, percent_tpm={}, experiment={}, source={}, date={})>".format(
            self.gene, self.tss, self.chrom, self.start, self.end, self.strand, self.cell_or_tissue, self.replicate, self.tpm, self.percent_tpm, self.experiment, self.source_name, self.date)