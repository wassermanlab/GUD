from sqlalchemy import (
    Column, Index, PrimaryKeyConstraint, String, ForeignKey,
    UniqueConstraint, CheckConstraint, Integer
)
from sqlalchemy.dialects import mysql
from GUD2.ORM.region import Region
from GUD2.ORM.source import Source
from GUD2.ORM.base import Base
from binning import containing_bins, contained_bins

class ShortTandemRepeat(Base):

    __tablename__ = "short_tandem_repeats"

    uid = Column("uid", mysql.INTEGER(unsigned=True))
    regionID = Column("regionID", Integer, ForeignKey('regions.uid'), nullable=False)
    sourceID = Column("sourceID", Integer, ForeignKey('sources.uid'), nullable=False)
    motif = Column("motif", String(30), nullable=False)
    pathogenicity = Column("pathogenicity", mysql.INTEGER(unsigned=True), nullable=False)

    __table_args__ = (
        PrimaryKeyConstraint(uid),
        UniqueConstraint(regionID, sourceID),

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
        # print(chrom, start, end)
        # q = Region.select_by_bin_range(session, chrom, start, end, [], True, False)

        bins = set(containing_bins(start, end) + contained_bins(start, end))
        
        q = session.query(cls, Region).\
        join().\
        filter(Region.uid == cls.regionID).\
        filter(Region.chrom == chrom, Region.end > start, Region.start < end).\
        filter(Region.bin.in_(bins))
        return q.all()
    
    @classmethod 
    def select_by_pathogenicity(cls, session):
        """returns all strs that are pathogenic"""
        q = session.query(cls).filter(cls.pathogenicity != 0)

        return q.all() 

    @classmethod
    def select_by_motif(cls, session, motif, compute_rotations=False):
        """returns all occurences of a certain motif computing 
        rotations if requested"""
        if compute_rotations is True: #compute the rotations
            motif_len = len(motif)
            motif_temp = motif 
            motifs = []
            for i in range(motif_len):
                motif_temp = motif_temp[1:motif_len] + motif_temp[0]
                motifs.append(motif_temp)    
        else:    
            motifs = [motif]

        q = session.query(cls).filter(cls.motif.in_(motifs))

        return q.all()    

    @classmethod
    def is_unique(cls, session, regionID, sourceID):
        q = session.query(cls).filter(cls.regionID == regionID, cls.sourceID == sourceID)
        q = q.all()
        if len(q) == 0:
            return True
        else: 
            return False 

    def __str__(self):
        return "{}\t{}".format(self.motif, self.pathogenicity)

    def __repr__(self):
        return "<ShortTandemRepeat(uid={}, regionID={}, sourceID={}, motif={}, pathogencity={})>".format(
            self.uid, self.regionID, self.sourceID, self.motif, self.pathogenicity)
