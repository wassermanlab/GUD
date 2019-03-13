from binning import (
    assign_bin
)
from sqlalchemy import (
    Column,
    ForeignKey,
    Index,
    Integer,
    PrimaryKeyConstraint,
    String,
    UniqueConstraint
)
from sqlalchemy.dialects import mysql

from .base import Base
from .region import Region
from .source import Source

class Gene(Base):

    __tablename__ = "genes"

    uid = Column(
        "uid",
        mysql.INTEGER(unsigned=True)
    )

    regionID = Column(
        "regionID",
        Integer,
        ForeignKey("regions.uid"),
        nullable=False
    )

    sourceID = Column(
        "sourceID",
        Integer,
        ForeignKey("sources.uid"),
        nullable=False
    )

    name = Column(
        "name",
        String(75),
        nullable=False
    )

    name2 = Column(
        "name2",
        String(75),
        nullable=False
    )

    cdsStart = Column(
        "cdsStart",
        mysql.INTEGER(unsigned=True),
        nullable=False
    )

    cdsEnd = Column(
        "cdsEnd",
        mysql.INTEGER(unsigned=True),
        nullable=False
    )

    exonStarts = Column(
        "exonStarts",
        mysql.LONGBLOB,
        nullable=False
    )

    exonEnds = Column(
        "exonEnds",
        mysql.LONGBLOB,
        nullable=False
    )

    __table_args__ = (
        PrimaryKeyConstraint(uid),
        UniqueConstraint(
            regionID,
            sourceID,
            name
        ),
        Index("ix_regionID", regionID), # query by bin range
        Index("ix_name", name),
        Index("ix_name2", name2),
        {
            "mysql_engine": "MyISAM",
            "mysql_charset": "utf8"
        }
    )

    @classmethod
    def is_unique(cls, session, regionID, sourceID,
        strand, name):

        q = session.query(cls).\
            filter(
                cls.regionID == regionID,
                cls.sourceID == sourceID,
                cls.strand == strand,
                cls.name == name
            )

        return len(q.all()) == 0

    @classmethod
    def select_unique(cls, session, regionID,
        sourceID, strand, name):

        q = session.query(cls).\
            filter(
                cls.regionID == regionID,
                cls.sourceID == sourceID,
                cls.strand == strand,
                cls.name == name
            )

        return q.first()

    @classmethod
    def select_by_location(cls, session, chrom,
        start, end):
        """
        Query objects by genomic location.
        """

        bin = assign_bin(start, end)

        q = session.query(cls, Region).join().\
            filter(Region.uid == cls.regionID).\
            filter(
                Region.chrom == chrom,
                Region.end > start,
                Region.start < end
            ).\
            filter(Region.bin == bin)

        return q.all()

    @classmethod
    def select_by_name(cls, session, name):
        """
        Query objects by gene symbol.
        """

        q = session.query(cls).\
            filter(
                cls.name2 == name
            )

        return q.all()

    @classmethod
    def select_by_names(cls, session, names=[]):
        """
        Query objects by multiple gene symbols.
        If no genes are provided, return all objects.
        """

        q = session.query(cls)

        if names:
            q = q.filter(cls.name2.in_(names))

        return q.all()

    @classmethod
    def select_by_uid(cls, session, uid):
        """
        Query objects by uid.
        """

        q = session.query(cls).\
            filter(cls.uid == uid)

        return q.first()

    @classmethod
    def select_by_uid_joined(cls, session, uid):
        """
        Query objects by uid.
        """

        q = session.query(cls, Region).join().\
            filter(Region.uid == cls.regionID).\
            filter(cls.uid == uid)

        return q.first()

    @classmethod
    def get_all_gene_symbols(cls, session):
        """
        Return the gene symbol (name2 field) of all
        objects.
        """

        q = session.query(cls.name2).distinct().all()

        return [g[0] for g in q]

    def __repr__(self):

        return "<Gene(%s, %s, %s, %s)>" % \
            (
                "uid={}".format(self.uid),
                "name={}".format(self.name),
                "name2={}".format(self.name2),
                "strand={}".format(self.strand)
            )