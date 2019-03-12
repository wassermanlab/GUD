from binning import (
    containing_bins,
    contained_bins
)
from sqlalchemy import (
    Column, Index, PrimaryKeyConstraint, String, ForeignKey,
    UniqueConstraint, CheckConstraint, Integer
)
from sqlalchemy.dialects import mysql

from .base import Base
from .region import Region
from .source import Source

class CNV(Base):
    __tablename__ = "copy_number_variants"

    uid = Column("uid", String(50))
    regionID = Column("regionID", Integer, ForeignKey('regions.uid'), nullable=False)
    sourceID = Column("sourceID", Integer, ForeignKey('sources.uid'), nullable=False)
    variant_type = Column("variant_type", String(50), nullable=False)
    copy_number = Column("copy_number", Integer, nullable=False)
    clinical_interpretation = Column("clinical_interpretation", String(50), nullable=False)

    __table_args__ = (
        PrimaryKeyConstraint(uid),

        Index("ix_cnv", regionID),
        Index("ix_cnv_uid", uid),
        Index("ix_cnv_clinical_interpretation", clinical_interpretation),

        {
            "mysql_engine": "MyISAM",
            "mysql_charset": "utf8"
        }
    )

    @classmethod
    def select_by_location(cls, session, chrom, start, end):
        """
        Query objects based off of their location being within the start only
        motifs through that  
         """
        bins = set(containing_bins(start, end) + contained_bins(start, end))
        q = session.query(cls, Region).\
        join().\
        filter(Region.uid == cls.regionID).\
        filter(Region.chrom == chrom, Region.end > start, Region.start < end).\
        filter(Region.bin.in_(bins))
        return q.all()

    @classmethod
    def select_by_exact_location(cls, session, chrom, start, end):
        """
        Query objects based off of their location being within the start only
        motifs through that  
         """
        bin = assign_bin(start, end)
        
        q = session.query(cls, Region).\
        join().\
        filter(Region.uid == cls.regionID).\
        filter(Region.chrom == chrom, Region.start == start, Region.end == end).\
        filter(Region.bin == bin)
        return q.first()  

    @classmethod
    def select_by_uid(cls, session, uid):
        q = session.query(cls, Region).\
        join().\
        filter(Region.uid == cls.regionID).\
        filter(cls.uid == uid)
        return q.first()

    @classmethod
    def is_unique(cls, session, name):
        q = session.query(cls).filter(cls.uid == name)
        q = q.all()
        return len(q) == 0

    def __str__(self):
        return "{}\t{}".format(self.copy_number, self.clinical_significance)

    def __repr__(self):
        return "<CNV(uid={}, regionID={}, sourceID={}, copy_number={}, clinical_significance={}, variant_type={}yes)>".format(
            self.uid, self.regionID, self.sourceID, self.copy_number, self.clinical_interpretation, self.variant_type)
