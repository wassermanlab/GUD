from lib.GUD.utils.bin_range import BinRange

from sqlalchemy import (
    Column, Index, PrimaryKeyConstraint, String
)

from sqlalchemy.dialects import mysql
from sqlalchemy.ext.declarative import declarative_base

Base = declarative_base()

class ChromSize(Base):
    
    __tablename__ = "chrom_size"

    chrom = Column("chrom", String(5), nullable=False)
    size = Column("size", mysql.INTEGER(unsigned=True), nullable=False)

    __table_args__ = (

        PrimaryKeyConstraint(
            chrom
        ),

        Index("ix_chrom_size", chrom),

        {
            "mysql_engine": "MyISAM",
            "mysql_charset": "utf8"
        }
    )

    @classmethod
    def chrom_sizes(cls, session): 
        """
        Returns the sizes of all chroms as a dict.
        """

        chrom_sizes = {}
        
        q = session.query(cls)
        
        for c in q:
            chrom_sizes.setdefault(c.chrom, int(c.size))

        return chrom_sizes

    @classmethod
    def chrom_size(cls, session, chrom):
        """
        Returns the size of the given chrom.
        """
        
        q = session.query(cls).filter(cls.chrom == chrom)

        for c in q:
            return int(c.size)

    def __str__(self):
        return "{}\t{}".format(self.chrom, self.size)

    def __repr__(self):
        return "<ChromSize(chrom={}, size={})>".format(
            self.chrom, self.size)
