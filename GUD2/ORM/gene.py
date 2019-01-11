import re

from sqlalchemy import (
    Column, Index, PrimaryKeyConstraint, String, ForeignKey,
    UniqueConstraint, CheckConstraint, Integer
)
from sqlalchemy.dialects import mysql
from GUD2.ORM.region import Region
from GUD2.ORM.source import Source
from GUD2.ORM.base import Base
from binning import containing_bins, contained_bins, assign_bin

class Gene(Base):

    __tablename__ = "genes"

    uid = Column("uid", mysql.INTEGER(unsigned=True))
    regionID = Column("regionID", Integer, ForeignKey('regions.uid'), nullable=False)
    sourceID = Column("sourceID", Integer, ForeignKey('sources.uid'), nullable=False)
    name = Column("name", String(75), nullable=False)
    cdsStart = Column("cdsStart", mysql.INTEGER(unsigned=True), nullable=False)
    cdsEnd = Column("cdsEnd", mysql.INTEGER(unsigned=True), nullable=False)
    exonStarts = Column("exonStarts", mysql.LONGBLOB, nullable=False)
    exonEnds = Column("exonEnds", mysql.LONGBLOB, nullable=False)
    name2 = Column("name2", String(75), nullable=False)
    strand = Column("strand", mysql.CHAR(1), nullable=False)

    __table_args__ = (
        PrimaryKeyConstraint(uid),
        UniqueConstraint(regionID, sourceID, strand, name),

        Index("ix_gene", regionID),
        Index("ix_gene_acce", name),
        Index("ix_gene_name", name2),

        {
            "mysql_engine": "MyISAM",
            "mysql_charset": "utf8"
        }
    )

    @classmethod
    def select_by_location(cls, session, chrom, start, end):
        """Query objects based off of their location being within the start only
        motifs through that  
        """
        bin = assign_bin(start, end)

        q = session.query(cls, Region).\
        join().\
        filter(Region.uid == cls.regionID).\
        filter(Region.chrom == chrom, Region.end > start, Region.start < end).\
        filter(Region.bin == bin)

        return q.all()

    @classmethod
    def select_by_name(cls, session, name):
        """Query refGene objects by common name. If no name provided,
        query all genes.
        """
        q = session.query(cls)

        if name:
            q = q.filter(cls.name2 == name)

        return q.all()

    @classmethod
    def select_by_names(cls, session, names=[]):
        """Query refGene objects by list of common names. If no names
        provided, query all genes.
        """
        q = session.query(cls)

        if names:
            q = q.filter(cls.name2.in_(names))

        return q.all()
    
    @classmethod
    def select_by_uid(cls, session, uid):
        """Query refGene objects by uid returning one uid"""
        q = session.query(cls)

        q.filter(cls.uid == uid)

        return q.first()
    
    @classmethod
    def select_by_uid_joined(cls, session, uid):
        """Query refGene objects by uid returning one uid"""
        q = session.query(cls, Region).\
        join().\
        filter(Region.uid == cls.regionID).\
        filter(cls.uid == uid)

        return q.first()

    @classmethod
    def get_all_gene_symbols(cls, session):
        """Query all gene objects in the database and return their gene
        symbols (name2 field).
        """
        q = session.query(cls.name2).distinct().all()

        return [g[0] for g in q]

    def __repr__(self):
        return "<Gene(uid={}, name={}, name2={}, strand={})>".format(
            self.uid, self.name, self.name2, self.strand,)

