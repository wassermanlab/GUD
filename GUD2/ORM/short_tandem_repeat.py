from sqlalchemy import (
    Column, Index, PrimaryKeyConstraint, String, ForeignKey,
    UniqueConstraint, CheckConstraint, Integer
)
from sqlalchemy.orm import relationship
from sqlalchemy.dialects import mysql
from GUD2.ORM.region import Region
from GUD2.ORM.source import Source
from GUD2.ORM.base import Base

class ShortTandemRepeat(Base):

    __tablename__ = "short_tandem_repeat"

    uid = Column("uid", mysql.INTEGER(unsigned=True))
    regionID = Column("regionID", Integer, ForeignKey('region.uid'), nullable=False)
    sourceID = Column("sourceID", Integer, ForeignKey('source.uid'), nullable=False)
    motif = Column("motif", String(30), nullable=False)
    pathogenicity = Column("pathogenicity", mysql.INTEGER(unsigned=True), nullable=False)

    __table_args__ = (
        PrimaryKeyConstraint(uid),
        UniqueConstraint(uid, regionID, sourceID),

        Index("ix_str", regionID),
        Index("ix_str_pathogenic", pathogenicity),
        Index("ix_str_motif", motif),

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
        region = Region()
        regions_in_range = region.select_by_bin_range(session, chrom, start, 
        end, [], True)

        q = session.query(cls, regions_in_range).\
        filter(cls.regionID == regions_in_range.uid)

        return q.all()

    def __str__(self):
        return "{}\t{}".format(self.motif, self.pathogenicit)

    def __repr__(self):
        return "<ShortTandemRepeat(uid={}, regionID={}, sourceID={}, motif={}, pathogencity={})>".format(
            self.uid, self.regionID, self.sourceID, self.motif, self.pathogenicity)
