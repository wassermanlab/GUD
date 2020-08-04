from sqlalchemy import (Column, Index, String, UniqueConstraint)
from sqlalchemy.dialects import mysql
from .genomicFeatureMixin1 import GFMixin1
from .base import Base
from .region import Region
from .source import Source
from sqlalchemy.ext.declarative import declared_attr


class Gene(GFMixin1, Base):
    # table declerations
    __tablename__ = "genes"

    # inherits uid, regionID, sourceID
    name = Column("name", String(75), nullable=False)
    gene_symbol = Column("gene_symbol", String(75), nullable=False)
    coding_start = Column("coding_start", mysql.INTEGER(unsigned=True), nullable=False)
    coding_end = Column("coding_end", mysql.INTEGER(unsigned=True), nullable=False)
    exon_starts = Column("exon_starts", mysql.LONGBLOB, nullable=False)
    exon_ends = Column("exon_ends", mysql.LONGBLOB, nullable=False)
    strand = Column("strand", mysql.CHAR(1), nullable=False)

    @declared_attr
    def __table_args__(cls):
        return (
            UniqueConstraint(cls.region_id, cls.name, cls.source_id, cls.strand),
            # query by bin range
            Index("ix_join", cls.source_id, cls.region_id),
            Index("ix_name", cls.name),
            Index("ix_gene_symbol", cls.gene_symbol),
            {"mysql_engine": "InnoDB", "mysql_charset": "utf8", }
        )

    # class methods
    @classmethod
    def select_by_names(cls, session, query, names=[]):
        """
        Query objects by multiple gene symbols.
        If no genes are provided, return all
        objects.
        """
        q = cls.make_query(session, query)
        q = q.filter(cls.gene_symbol.in_(names))
        return q

    @classmethod
    def get_all_gene_symbols(cls, session):
        """
        Return the gene symbol of all objects.
        """
        q = session.query(cls.gene_symbol).distinct()
        return q

    @classmethod
    def is_unique(cls, session, regionID, sourceID, name, strand):
        """
        Checks uniqueness by region source and name
        """
        q = session.query(cls).filter(cls.region_id == regionID, cls.strand == strand, 
                                      cls.name == name, cls.source_id == sourceID,)
        return len(q.all()) == 0

    @classmethod
    def as_genomic_feature(self, feat):
        """
        extend parent class by adding qualifiers
        """
        # parse the exon starts and ends
        exonStarts = []
        exonEnds = []

        # For each exon start...
        for i in str(feat.Gene.exon_starts).split(","):
            if i.isdigit():
                exonStarts.append(int(i))

        # For each exon end...
        for i in str(feat.Gene.exon_ends).split(","):
            if i.isdigit():
                exonStarts.append(int(i))

        # Define qualifiers
        qualifiers = {
            "uid": feat.Gene.uid,
            "name": feat.Gene.name,
            "gene_symbol": feat.Gene.gene_symbol,
            "coding_start": int(feat.Gene.coding_start),
            "coding_end": int(feat.Gene.coding_end),
            "exon_starts": feat.Gene.exon_starts,
            "exon_ends": feat.Gene.exon_ends,
            "source": feat.sourceName,
        }
        genomic_feature = super().as_genomic_feature(feat)
        genomic_feature.qualifiers = qualifiers
        return genomic_feature