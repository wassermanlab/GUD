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
    name2 = Column("name2", String(75), nullable=False)
    cdsStart = Column("cdsStart", mysql.INTEGER(unsigned=True), nullable=False)
    cdsEnd = Column("cdsEnd", mysql.INTEGER(unsigned=True), nullable=False)
    exonStarts = Column("exonStarts", mysql.LONGBLOB, nullable=False)
    exonEnds = Column("exonEnds", mysql.LONGBLOB, nullable=False)
    strand = Column("strand", mysql.CHAR(1))

    @declared_attr
    def __table_args__(cls):
        return (
            UniqueConstraint(cls.region_id, cls.name, cls.source_id),
            # query by bin range
            Index("ix_source_join", cls.source_id, cls.region_id),
            Index("ix_name", cls.name),
            Index("ix_name2", cls.name2),
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
        q = q.filter(cls.name2.in_(names))
        return q

    @classmethod
    def get_all_gene_symbols(cls, session):
        """
        Return the gene symbol (name2 field) of all objects.
        """
        q = session.query(cls.name2).distinct()
        return q

    @classmethod
    def is_unique(cls, session, regionID, sourceID, name):
        """
        Checks uniqueness by region source and name
        """
        q = session.query(cls).filter(cls.region_id == regionID,
                                      cls.name == name, cls.source_id == sourceID)
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
        for i in str(feat.Gene.exonStarts).split(","):
            if i.isdigit():
                exonStarts.append(int(i))

        # For each exon end...
        for i in str(feat.Gene.exonEnds).split(","):
            if i.isdigit():
                exonStarts.append(int(i))

        # Define qualifiers
        qualifiers = {
            "uid": feat.Gene.uid,
            "accession_number": feat.Gene.name,
            "gene_symbol": feat.Gene.name2,
            "cdsStart": int(feat.Gene.cdsStart),
            "cdsEnd": int(feat.Gene.cdsEnd),
            "exonStarts": feat.Gene.exonStarts,
            "exonEnds": feat.Gene.exonEnds,
            "source": feat.Source.name,
        }
        genomic_feature = super().as_genomic_feature(feat)
        genomic_feature.qualifiers = qualifiers
        return genomic_feature
