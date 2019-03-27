from sqlalchemy import (
    Column,
    Index,
    PrimaryKeyConstraint,
    String
)
from sqlalchemy.dialects import mysql

from .base import Base

class Chrom(Base):

    __tablename__ = "chroms"

    chrom = Column(
        "chrom",
        String(5),
        nullable=False
    )

    size = Column(
        "size",
        mysql.INTEGER(unsigned=True),
        nullable=False
    )

    __table_args__ = (
        PrimaryKeyConstraint(chrom),
        Index("ix_chrom", chrom),
        {
            "mysql_engine": "MyISAM",
            "mysql_charset": "utf8"
        }
    )

    @classmethod
    def select_by_chrom(cls, session, chrom):
        """
        Query objects by chromosome name.
        """

        q = session.query(cls)\
            .filter(cls.chrom == chrom)

        return q.first()

    @classmethod
    def select_by_chroms(cls, session, chroms=[]):
        """
        Query objects by multiple chromosome names.
        If no names are provided, return all
        objects.
        """

        q = session.query(cls)

        if chroms:
            q = q.filter(cls.chrom.in_(chroms))

        return q.all()

    @classmethod
    def chrom_size(cls, session, chrom):
        """
        Return the size of the given chrom.
        """
        
        q = cls.select_by_chroms(session, chrom)

        for c in q:
            return int(c.size)

    @classmethod
    def chrom_sizes(cls, session, chroms=[]): 
        """
        Return the size of the given chroms as a
        dict. If no chroms are provided, return
        the size of all chroms.
        """

        sizes = {}

        q = cls.select_by_chroms(session, chroms)
        
        for c in q:
            sizes.setdefault(c.chrom, int(c.size))

        return sizes

    def __repr__(self):

        return "<Chrom(%s, %s)>" % \
            (
                "chrom={}".format(self.chrom),
                "size={}".format(self.size)
            )