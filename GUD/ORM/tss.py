from GUD.bin_range import BinRange

from Bio.SeqFeature import FeatureLocation

from sqlalchemy import (
    Table, MetaData, Column, Integer, Float, CHAR, String, Enum, BLOB, PrimaryKeyConstraint, UniqueConstraint, Index, ForeignKey, join,and_
)

from sqlalchemy.dialects import mysql
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship, column_property

Base = declarative_base()

# XXX we should define this in a better way (i.e. in the database itself)
__GENE_ID_TYPE_HGNC_SYMBOL__ = 4

class Expression(Base):
    """ This is modelled after the CAGEd-oPOSSUM 'expression' table.

    This follows the SQLAlchemy model of an 'association' table defining a
    many-to-many relationship between fantom5_cage_peaks and fantom5_samples.
    """

    __tablename__ = 'fantom5_expression'

    cage_peak_id = Column('cage_peak_id', mysql.INTEGER(unsigned=True),
                          ForeignKey("fantom5_cage_peaks.id"),
                          primary_key=True)
    sample_id = Column('sample_id', mysql.INTEGER(unsigned=True),
                       ForeignKey("fantom5_samples.id"),
                       primary_key=True)
    tag_count = Column('tag_count', mysql.INTEGER(unsigned=True),
                       nullable=False)
    tpm = Column('tpm', mysql.DOUBLE(unsigned=True), nullable=False)
    relative_expression = Column('relative_expression',
                                  mysql.DOUBLE(unsigned=True), nullable=False)

    #sample = relationship("Sample", back_populates="fantom5_expression")
    #cage_peak = relationship("CagePeak", back_populates="fantom5_expression")
    sample = relationship("Sample", backref="fantom5_expression")
    cage_peak = relationship("CagePeak", backref="fantom5_expression")


    __table_args__ = (
        Index(sample_id),

        {
            'mysql_engine':'MyISAM',
            'mysql_charset':'utf8'
        }
    )


    def __str__(self):
        #return "{}\t{}\t{}\t{}\t{}".format(self.cage_peak_id, self.sample_id, self.tag_count, self.tpm, self.relative_expression)
        return "{}\t{}\t{}".format(self.tag_count, self.tpm, self.relative_expression)

    def __repr__(self):
        return "<Expression(cage_peak_id={}, sample_id={}, tag_count={}, tpm={}, relative_expression={})>".format(self.cage_peak_id, self.sample_id, self.tag_count, self.tpm, self.relative_expression)


class Sample(Base):
    """ This is modelled after the CAGEd-oPOSSUM 'experiments' table.
    """

    __tablename__ = 'fantom5_samples'

    id = Column('id', mysql.INTEGER(unsigned=True), nullable=False)
    FF_id = Column('FF_id', CHAR(12), nullable=False)
    CNhs_id = Column('CNhs_id', mysql.INTEGER(unsigned=True), nullable=False)
    type = Column('type', String(20))
    method = Column('method', String(12))
    name = Column('name', String(255), nullable=False)

    #cage_peaks = relationship("Expression", backref="famtom5_samples")
    expression = relationship("Expression", backref="famtom5_samples")

    #cage_peaks = relationship("CagePeak",
    #                          secondary="fantom5_expression",
    #                          back_populates="famtom5_samples")


    __table_args__ = (
        PrimaryKeyConstraint(id),
        Index(FF_id),

        {
            'mysql_engine':'MyISAM',
            'mysql_charset':'utf8'
        }
    )


    def select_cage_peak_expression(self, session, tag_count=0, tpm=0.0,
                                    relative_expression=0.0):
        """ Query all CAGE peaks that have the given level of expression for
        this sample.
        """

        q = (session.query(CagePeak, Expression)
                .join(Expression)
                .filter(Expression.sample_id == self.id,
                        CagePeak.id == Expression.cage_peak_id)
            )

        if tag_count:
            q = q.filter(Expression.tag_count >= tag_count)

        if tpm:
            q = q.filter(Expression.tpm >= tpm)

        if relative_expression:
            q = q.filter(Expression.relative_expression >= relative_expression)

        return q.all()


    def __str__(self):
        return "FF:{}\tCNhs{}\t{}\t{}\t{}".format(self.FF_id, self.CNhs_id, self.type, self.method, self.name)

    def __repr__(self):
        return "<Sample(id={}, FF_id={}, CNhs_id={}, type={}, method={}, name={})>".format(self.id, self.FF_id, self.CNhs_id, self.type, self.method, self.name)


class CagePeak(Base):
    """ This is actually modelled after the CAGEd-oPOSSUM 'tss' table.
    """

    __tablename__ = 'fantom5_cage_peaks'

    id = Column('id', mysql.INTEGER(unsigned=True), nullable=False)
    #
    # In the corresponding CAGEd-oPOSSUM table there is a search_region_id
    # which is essentially a bin number which is dynamically computed based
    # on merging CAGE peak regions with fixed length flanks into larger
    # regions. We will replace this with a bin field computed in the standard
    # UCSC fashion.
    #
    bin = Column('bin', mysql.SMALLINT(unsigned=True), nullable=False)
    name = Column('name', CHAR(30), nullable=False)
    chrom = Column('chrom', CHAR(5), nullable=False)
    start = Column('start', mysql.INTEGER(unsigned=True), nullable=False)
    end = Column('end', mysql.INTEGER(unsigned=True), nullable=False)
    strand = Column('strand', CHAR(1))
    is_tss = Column('is_tss', mysql.TINYINT(), nullable=False)
    max_tag_count = Column('max_tag_count', mysql.INTEGER(unsigned=True),
                            nullable=False)
    max_tpm = Column('max_tpm', mysql.DOUBLE())

    #samples = relationship("Expression", backref="fantom5_cage_peaks")
    expression = relationship("Expression", backref="famtom5_cage_peaks")

    #samples = relationship("Sample",
    #                       secondary="fantom5_expression",
    #                       back_populates="famtom5_cage_peaks")

    #cage_peak_genes = relationship("CagePeakGene",
    #                               back_populates="fantom5_cage_peaks")


    __table_args__ = (
        PrimaryKeyConstraint(id),
        UniqueConstraint(name),
        Index(chrom, bin, start),

        {
            'mysql_engine':'MyISAM',
            'mysql_charset':'utf8'
        }
    )

    def select_sample_expression(self, session, tag_count=0, tpm=0.0,
                                 relative_expression=0.0):
        """ Query all samples that have the given level of expression for
        this cage peak.
        """

        q = (session.query(Sample, Expression)
                .join(Expression)
                .filter(Expression.cage_peak_id == self.id,
                        Sample.id == Expression.sample_id)
            )

        if tag_count:
            q = q.filter(Expression.tag_count >= tag_count)

        if tpm:
            q = q.filter(Expression.tpm >= tpm)

        if relative_expression:
            q = q.filter(Expression.relative_expression >= relative_expression)

        return q.all()

    @classmethod
    def select_by_bin_range(cls, session, chrom, start, end, strand=None,
                            is_tss=False, bins=[], compute_bins=True): 
        """ Query objects by chromosomal range using the binning system to
        speed up range searches. If bins are provided use the given bins.
        If bins are NOT provided AND compute_bins is set to True, then
        compute the bins. Otherwise perform the range query without the use
        of bins.
        """

        if not bins and compute_bins:
            br = BinRange()
            bins = br.allBinsInRange(start, end)

        q = session.query(cls).filter(
                cls.chrom == chrom, cls.end > start, cls.start < end)

        if strand:
            if strand == 1:
                strand == '+'
            elif strand == -1:
                strand == '-'

            if strand != '+' and strand != '-':
                raise ValueError("Strand should be 1, -1, '+', '-' or None")

            q = q.filter(cls.strand == strand)

        if bins:
            q = q.filter(cls.bin.in_(bins))

        return q.all()

    def __str__(self):
        return "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(self.name, self.chrom, self.start, self.end, self.strand, self.is_tss, self.max_tag_count, self.max_tpm)

    def __repr__(self):
        return "<CagePeak(id={}, name={}, chrom={}, start={}, end={}, strand={}, is_tss={}, max_tag_count={}, max_tpm={})>".format(self.id, self.name, self.chrom, self.start, self.end, self.strand, self.is_tss, self.max_tag_count, self.max_tpm)

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            # Should we just use ID, just name, or combination of
            # chrom/start/end/strand, or all of the above?
            return self.name == other.name
        return NotImplemented

    def __ne__(self, other):
        if isinstance(other, self.__class__):
            return not self.__eq__(other)
        return NotImplemented


class CagePeakGene(Base):
    """
    This table describes the association of genes to CAGE peaks.
    """

    __tablename__ = 'fantom5_cage_peak_genes'

    cage_peak_id = Column('cage_peak_id', mysql.INTEGER(unsigned=True),
                          ForeignKey("fantom5_cage_peaks.id"),
                          primary_key=True)
    gene_id_type = Column('gene_id_type', mysql.TINYINT(unsigned=True),
                          primary_key=True)
    gene_id = Column('gene_id', CHAR(16), primary_key=True)

    #cage_peak = relationship("CagePeak",
    #                         back_populates="fantom5_cage_peak_genes")


    __table_args__ = (
        Index(gene_id),

        {
            'mysql_engine':'MyISAM',
            'mysql_charset':'utf8'
        }
    )


    def __str__(self):
        return "{}\t{}\t{}".format(self.cage_peak_id, self.gene_id_type, self.gene_id)

    def __repr__(self):
        return "<CagePeakGene(cage_peak_id={}, gene_id_type={}, gene_id={})>".format(self.cage_peak_id, self.gene_id_type, self.gene_id)


class Enhancer(Base):
    """ FANTOM5 enhancer object. Similar to a CagePeak.
    """

    __tablename__ = 'fantom5_enhancers'

    #id = Column('id', mysql.INTEGER(unsigned=True), nullable=False)
    bin = Column('bin', mysql.SMALLINT(unsigned=True), nullable=False)
    #name = Column('name', CHAR(30), nullable=False)
    chrom = Column('chrom', CHAR(5), nullable=False)
    start = Column('start', mysql.INTEGER(unsigned=True), nullable=False)
    end = Column('end', mysql.INTEGER(unsigned=True), nullable=False)
    cell_or_tissue = Column('cell_or_tissue', String(40), nullable=False)
    #score = Column('score', mysql.INTEGER(unsigned=True))


    __table_args__ = (
        PrimaryKeyConstraint(chrom, cell_or_tissue, start, end),
        #UniqueConstraint(name),
        Index(chrom, bin, start),

        {
            'mysql_engine':'MyISAM',
            'mysql_charset':'utf8'
        }
    )

    @classmethod
    def select_by_bin_range(cls, session, chrom, start, end, cell_or_tissue=[], strand=None,
                            bins=[], compute_bins=True): 
        """ Query objects by chromosomal range using the binning system to
        speed up range searches. If bins are provided use the given bins.
        If bins are NOT provided AND compute_bins is set to True, then
        compute the bins. Otherwise perform the range query without the use
        of bins.
        """

        if not bins and compute_bins:
            br = BinRange()
            bins = br.allBinsInRange(start, end)

        q = session.query(cls).filter(
                cls.chrom == chrom, cls.end > start, cls.start < end)

        if bins:
            q = q.filter(cls.bin.in_(bins))

        if cell_or_tissue:
            q = q.filter(cls.cell_or_tissue.in_(cell_or_tissue))

        return q.all()

    def __str__(self):
        return "{}\t{}\t{}\t{}".format(self.chrom, self.start, self.end, self.cell_or_tissue)

    def __repr__(self):
        return "<Enhancer(chrom={}, start={}, end={}, cell={})>".format(self.chrom, self.start, self.end, self.cell_or_tissue)

#    def __eq__(self, other):
#        if isinstance(other, self.__class__):
#            # Should we just use ID, just name, or combination of
#            # chrom/start/end/strand, or all of the above?
#            return self.name == other.name
#        return NotImplemented

#    def __ne__(self, other):
#        if isinstance(other, self.__class__):
#            return not self.__eq__(other)
#        return NotImplemented


class Fantom5(object):
    """ This is an attempt to encapsulate the 3-way join of CagePeak,
    Expression and Sample into a single FANTOM5 class.
    """

#    samp_cp_expr_join = (
#        join(Sample, Expression).join(Expression, CagePeak)
#    )
#
#    __table__ = samp_cp_expr_join
#
#    bin = column_property(CagePeak.bin)
#    FF_id = column_property(Sample.FF_id)
#    CNhs_id = column_property(Sample.CNhs_id)
#    type = column_property(Sample.type)
#    method = column_property(Sample.method)
#    sample_name = column_property(Sample.name)
#
#    tag_count = column_property(Expression.tag_count)
#    tpm = column_property(Expression.tpm)
#    relative_expression = column_property(Expression.relative_expression)
#
#    cage_peak_name_name = column_property(CagePeak.name)
#    chrom = column_property(CagePeak.chrom)
#    start = column_property(CagePeak.start)
#    end = column_property(CagePeak.end)
#    strand = column_property(CagePeak.strand)
#    is_tss = column_property(CagePeak.is_tss)
#    max_tag_count = column_property(CagePeak.max_tag_count)
#    max_tpm = column_property(CagePeak.max_tpm)


    @classmethod
    def select_by_bin_range(
            cls, session, chrom, start, end, ff_ids=[], tss_only=0,
            tag_count=0, tpm=0.0, relative_expression=0.0, bins=[],
            compute_bins=False): 

#        q = (session.query(cls)
#                .filter(cls.chrom == chrom,
#                        cls.end > start,
#                        cls.start < end)
#            )
#
#        if not bins and compute_bins:
#            br = BinRange()
#            bins = br.allBinsInRange(start, end)
#
#        if bins:
#            q = q.filter(cls.bin.in_(bins))
#
#        if ff_ids:
#            q = q.filter(cls.FF_id.in_(ff_ids))
#
#        if tss_only:
#            q = q.filter(cls.is_tss == 1)
#
#        if tag_count:
#            q = q.filter(cls.tag_count >= tag_count)
#
#        if tpm:
#            q = q.filter(cls.tpm >= tpm)
#
#        if relative_expression:
#            q = q.filter(cls.relative_expression >= relative_expression)
#
#        return q.all()

        q = (session.query(Expression, CagePeak, Sample)
                .join(CagePeak)
                .join(Sample)
                .filter(CagePeak.chrom == chrom,
                        CagePeak.end > start,
                        CagePeak.start < end)
            )

        if not bins and compute_bins:
            br = BinRange()
            bins = br.allBinsInRange(start, end)

        if bins:
            q = q.filter(CagePeak.bin.in_(bins))

        if ff_ids:
            q = q.filter(Sample.FF_id.in_(ff_ids))

        if tss_only:
            q = q.filter(CagePeak.is_tss == 1)

        if tag_count:
            q = q.filter(Expression.tag_count >= tag_count)

        if tpm:
            q = q.filter(Expression.tpm >= tpm)

        if relative_expression:
            q = q.filter(Expression.relative_expression >= relative_expression)

        feature_tuple = q.all()

        results = []
        for ft in feature_tuple:
            (expr, cp, samp) = ft

            results.append(
                Fantom5(samp.FF_id, samp.name, cp.name, cp.chrom,
                cp.start, cp.end, cp.strand, expr.tag_count, expr.tpm,
                expr.relative_expression)
            )

        return results


    def __init__(self, FF_id, sample_name, cage_peak_name, chrom, start, end, strand, tag_count, tpm, relative_expression):
        self.FF_id = FF_id
        self.sample_name = sample_name
        self.cage_peak_name = cage_peak_name
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.tag_count = tag_count
        self.tpm = tpm
        self.relative_expression = relative_expression

    def __str__(self):
        return "FF:{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.3f}\t{:.3f}".format(self.FF_id, self.sample_name, self.cage_peak_name, self.chrom, self.start, self.end, self.strand, self.tag_count, self.tpm, self.relative_expression)

    def __repr__(self):
        return "<Fantom5({}, {}, {}, {}, {}, {}, {}, {}, {:.3f}, {:.3f})>".format(self.FF_id, self.sample_name, self.cage_peak_name, self.chrom, self.start, self.end, self.strand, self.tag_count, self.tpm, self.relative_expression)


def select_by_bin_range(session, chrom, start, end, ff_ids=[],
        tss_only=0, tag_count=0, tpm=0.0, relative_expression=0.0,
        bins=[], compute_bins=False): 
    """ Module level query to fetch information joined across the three
    FANTOM5 tables/ORM classes allowing filtering with optional criteria
    on each of the 3 tables. If bins are provided use the given bins.
    If bins are NOT provided AND compute_bins is set to True, then
    compute the bins. Otherwise perform the range query without the use
    of bins.
    """

    q = (session.query(Expression, CagePeak, Sample)
            .join(CagePeak)
            .join(Sample)
            .filter(CagePeak.chrom == chrom,
                    CagePeak.end > start,
                    CagePeak.start < end)
        )

    if not bins and compute_bins:
        br = BinRange()
        bins = br.allBinsInRange(start, end)

    if bins:
        q = q.filter(CagePeak.bin.in_(bins))

    if ff_ids:
        q = q.filter(Sample.FF_id.in_(ff_ids))

    if tss_only:
        q = q.filter(CagePeak.is_tss == 1)

    if tag_count:
        q = q.filter(Expression.tag_count >= tag_count)

    if tpm:
        q = q.filter(Expression.tpm >= tpm)

    if relative_expression:
        q = q.filter(Expression.relative_expression >= relative_expression)

    feature_tuple = q.all()

    results = []
    for ft in feature_tuple:
        (expr, cp, samp) = ft

        results.append(
            Fantom5(samp.FF_id, samp.name, cp.name, cp.chrom,
            cp.start, cp.end, cp.strand, expr.tag_count, expr.tpm,
            expr.relative_expression)
        )

    return results


def select_ref_gene_cage_peaks(session, ref_gene, by_association=False,
        by_proximity=False, same_strand=True, is_tss=False, upstream_bp=0): 

    """ Module level query to query CAGE peaks by associated refSeq gene.

    The CAGE peaks are selected either by one of two methods:
    1) By association - all CAGE peaks which are explicitly identified as
       being associated to the gene by FANTOM5 are returned;
    OR
    2) By proximity - all CAGE peaks in proximity to the gene are returned.
       The region extends from the 3' end of the gene to the 5' end plus
       whatever amount of upstream bp are provided. Only CAGE peaks with the
       same strand are returned unless same_strand is set to False.
    """

    if by_proximity is False and by_association is False:
        raise ValueError(
            "Neither by_proximity nor by_association was specified")

    if by_proximity is True and by_association is True:
        raise ValueError(
            "Please specify either by_proximity or by_association but not both")

    if by_association is True:
        q = (session.query(CagePeak)
                .join(CagePeakGene)
                #.join(Expression)
            .filter(
            CagePeakGene.gene_id_type == __GENE_ID_TYPE_HGNC_SYMBOL__,
            CagePeakGene.gene_id == ref_gene.name2))

        if is_tss:
            q = q.filter(CagePeak.is_tss == 1)

        #if (tag_count):
        #    q = q.filter(Expression.tag_count >= tag_count)
        #if (tpm):
        #    q = q.filter(Expression.tpm >= tpm)
        #if (relative_expression):
        #    q = q.filter(Expression.relative_expression >= relative_expression)

        return q.all()

    elif by_proximity is True:
        strand = None
        if same_strand:
            strand = ref_gene.strand

        start = ref_gene.txStart
        end = ref_gene.txEnd
        if upstream_bp and upstream_bp > 0:
            if strand == 1 or strand == '+':
                start -= upstream_bp
            elif strand == -1 or strand == '-':
                end += upstream_bp

        return CagePeak.select_by_bin_range(
                    session, ref_gene.chrom, start, end, strand,
                    is_tss=is_tss, bins=None, compute_bins=True)



def select_sample_by_ff_id(session, ff_id):
    """ Fetch a sample by Fantom5 ID
    """

    q = session.query(Sample).filter(Sample.FF_id == ff_id)

    return q.one_or_none()


def get_ff_id(session, tissue):
    """ Fetch the FF_id from the given sample name (will search for wildcards on both sides)
    """

    #in order to use the LIKE command in SQL, need to append '%' as a wildcard flag
    #(% acts as * in shell script pattern matching)
    #tissue = "%"+tissue+"%"

    q = session.query(Sample.FF_id).filter(Sample.name.ilike(tissue))

    #return a list (of tuples -- format [(id#1,),(id#2,)]). May be empty, contain 1 element, or multiple elements
    return q.all()

def tss_by_ff_id(session, ffid, min_tpm):
    """Get the TSS ID and min tpm to get CAGE peak (TSSs) from that sample with at least the minimum tpm as specified
    """

    #assume that ffids are a list
    q = session.query(Sample.id).filter(Sample.FF_id.in_(ffid))
    #return list of all sample ids
    tuple_ids = q.all()
    #post-processing on this list, as all() returns a list of tuples in the form [(id#,),(id#,)..]
    ids = []
    for term in tuple_ids:
        #need to force type<int> as values from table will return long type with a trailing L
        ids.append(int(term[0]))
    #next query: filter Expression table on all sample_ids and min_tpm
    
#    q = session.query(Expression.cage_peak_id, Expression.tpm).filter(Expression.sample_id.in_(ids), Expression.tpm >= min_tpm)
    q = session.query(Expression.cage_peak_id, Expression.tpm).filter(Expression.sample_id.in_(ids))
    
    q = q.filter(Expression.tpm >= min_tpm)     
##
##    #since sample --> cage peak (tss) requires the sample id (not ff_id), need to do this query first
##    q = session.query(Sample.id).filter(Sample.FF_id == ffid)
##
##    # only expect a single return value, so the use of scalar() forces the method one() but returns a 'scalar' value rather than a list
##    sample_id = q.scalar()
##
##    #replacing q: as a standard, the variable q is used for all queries
##    #filtering on both sample and minimum tpm
##    q = session.query(Expression.cage_peak_id, Expression.tpm).filter(Expression.sample_id == sample_id, Expression.tpm >= min_tpm)

    #possible for one sample to have multiple assigned TSSs, so we will have to return all
    #output is list of tuples containing [(cage peak id/tss id, tpm),..] OR empty list
    return q.all()

def cage_peak_to_gene(session, cage_peak):
    """ Fetch gene name associated to the CAGE peak """

    q = session.query(CagePeakGene.gene_id).filter(CagePeakGene.gene_id_type == __GENE_ID_TYPE_HGNC_SYMBOL__, CagePeakGene.cage_peak_id == cage_peak)

    #each cage peak should be associated to only one gene, but there are some overlapping genes in F5 that associate the cage peak to
    #more than one gene. In this case, will need to fetch all, even though it will mostly be one (or no) associated RefSeq genes
    result = q.all() # returns list in form [(gene_name,)]
    #in-place method to return a list containing a single [gene], a list of [None], or a list of multiple [gene1,gene2]
    gene = []
    for item in result:
        gene.append(item[0])
    return gene
    #return q.one_or_none()

#def __str__(feature_tuple):
#    """ XXX Expects the input argument to be a tuple of a FANTOM5 Expression,
#    CagePeak and Sample objects. This DOES NOT REALLY WORK in the way we would
#    desire. What is desired is a to be able to call print on the fantom5
#    module level and have the this string method automatically invoked with
#    a python build in print function. This may require making a FANTOM5 class
#    within this module.
#    XXX
#    """
#    (expr, cp, samp) = feature_tuple
#
#    return "FF:{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.3f}\t{:.3f}".format(samp.FF_id, samp.name, cp.name, cp.chrom, cp.start, cp.end, cp.strand, expr.tag_count, expr.tpm, expr.relative_expression)
