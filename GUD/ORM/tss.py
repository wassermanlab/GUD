import re

from binning import containing_bins, contained_bins 
from Bio.SeqFeature import FeatureLocation

from sqlalchemy import (
    and_, or_, Column, Date, Enum, Float, Index,
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
    avg_tpm = Column("avg_tpm", Float, nullable=False)
#    percent_tpm = Column("percent_tpm", Float, nullable=False)
    experiment_type = Column("experiment_type", String(25), nullable=False)
    source_name = Column("source_name", String(25), nullable=False)
    date = Column("date", Date(), nullable=True)

    __table_args__ = (

        PrimaryKeyConstraint(
            chrom, start, end, strand, cell_or_tissue,
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
    def select_by_gene(cls, session, gene):
        """
        Query objects by gene. If no gene is provided, query all genes.
        """

        q = session.query(cls)

        if gene:
            q = q.filter(cls.gene == gene)

        return q.all()

    @classmethod
    def select_by_genes(cls, session, genes=[]):
        """
        Query objects by list of genes. If no genes are provided, query
        all genes.
        """

        q = session.query(cls)

        if genes:
            q = q.filter(cls.gene.in_(genes))

        return q.all()

    @classmethod
    def select_by_tss(cls, session, tss=[]):
        """
        Query objects by TSS. If no TSS is provided, query all TSSs.
        TSS is provided as a {list} of {lists}/{tuples} in the form 
        gene, tss.
        """

        q = session.query(cls)

        if gene and tss:
            q = q.filter(cls.gene == gene, cls.tss == tss)

        return q.all()

    @classmethod
    def select_by_multiple_tss(cls, session, tss=[]):
        """
        Query objects by list of TSSs. If no TSS are provided, query
        all TSSs. TSSs are provided as a {list} of {lists}/{tuples}
        in the form gene, tss.
        """

        q = session.query(cls)

        if tss:
            # Initialize
            ands = []
            # For each gene, TSS pair...
            for i, j in tss:
                ands.append(and_(cls.gene == i,
                    cls.tss == j))
            q = q.filter(or_(*ands))

        return q.all()

    @classmethod
    def select_by_sample(cls, session, sample=[]):
        """
        Query objects by list of samples. If no samples are provided,
        query all TSSs.
        """

        q = session.query(cls)
        
        if sample:
            q = q.filter(cls.cell_or_tissue.in_(sample))

        return q.all()

    @classmethod
    def select_by_differential_expression(cls, session, sample=[],
        avg_tpm=10.0, perc_tpm=25.0, expression_in_all_samples=False):
        """
        Query objects differentially expressed in sample.
        """

        # Initialize
        expression = {}
        float_regexp = re.compile("\d+\.\d+")

        # For each feat...
        for feat in cls.select_by_sample(session, sample=sample):
            # If not expression for TSS...
            if (feat.gene, feat.tss) not in expression:
                expression.setdefault((feat.gene, feat.tss), {})
                # For each sample...
                for sample in samples:
                    expression[(feat.gene, feat.tss)].setdefault(sample, 0.0)
            # Get TPMs
            tpms = map(float, re.findall(float_regexp, feat.tpm))
            expression[(feat.gene, feat.tss)][feat.cell_or_tissue] = sum(tpms) / len(tpms)
            # If enough TPMs...
            if sum(tpms) / len(tpms) > avg_tpm:
                tss.append((feat.gene, feat.tss))

        print(expression)
        exit(0)
#                
#        q = q.query(cls.cell_or_tissue.in_(sample))
#
#        q = session.query(cls, func.avg(cls.tpm).label("avg_tpm"),
#            func.sum(cls.percent_tpm).label("sum_perc_tpm")).group_by(
#            cls.gene, cls.tss, cls.chrom, cls.start, cls.end, cls.strand).having(
#            func.avg(cls.tpm) >= avg_tpm_thresh).having(
#            func.sum(cls.percent_tpm) >= sum_perc_tpm_thresh)
#
#        if sample:
#            q = q.filter(cls.cell_or_tissue.in_(sample))
#
#        q = q
#        print(q.all())
#        exit(0)

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
        return "{}\t{}\t{}\t{} ({})\t{}\t{}\t{}".format(self.chrom,
            self.start, self.end, self.gene, self.tss, self.strand,
            self.tpm, self.cell_or_tissue)

    def __repr__(self):
        return "<TSS(gene={}, tss={}, chrom={}, start={}, end={}, strand={}, sample={}, tpm={}, experiment={}, source={})>".format(
            self.gene, self.tss, self.chrom, self.start, self.end, self.strand, 
            self.cell_or_tissue, self.replicate, self.tpm, self.experiment_type,
            self.source_name)