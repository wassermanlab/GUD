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
from .genomic_feature import GenomicFeature

class ShortTandemRepeat(Base):

    __tablename__ = "short_tandem_repeats"

    uid = Column("uid", mysql.INTEGER(unsigned=True))
    regionID = Column("regionID", Integer, ForeignKey('regions.uid'), nullable=False)
    sourceID = Column("sourceID", Integer, ForeignKey('sources.uid'), nullable=False)
    motif = Column("motif", String(30), nullable=False)
    pathogenicity = Column("pathogenicity", mysql.INTEGER(unsigned=True), nullable=False)

    __table_args__ = (
        PrimaryKeyConstraint(uid),
        UniqueConstraint(regionID, sourceID, pathogenicity),

        Index("ix_str", regionID),
        Index("ix_str_pathogenic", pathogenicity),
        Index("ix_str_motif", motif),

        {
            "mysql_engine": "MyISAM",
            "mysql_charset": "utf8"
        }
    )

    @classmethod
    def select_by_location(cls, session, chrom, start, end, as_genomic_feature=False):
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
        if as_genomic_feature:
            feats = []
            # For each feature...
            for feat in q.all():
                feats.append(
                    cls.__as_genomic_feature(feat)
                )
            return feats
        return q.all()

    @classmethod
    def select_by_exact_location(cls, session, chrom, start, end, as_genomic_feature=False):
        """
        Query objects based off of their location being within the start only
        motifs through that  
         """
        # print(chrom, start, end)
        # q = Region.select_by_bin_range(session, chrom, start, end, [], True, False)

        bin = assign_bin(start, end)
        
        q = session.query(cls, Region).\
        join().\
        filter(Region.uid == cls.regionID).\
        filter(Region.chrom == chrom, Region.start == start, Region.end == end).\
        filter(Region.bin == bin)
        if as_genomic_feature:
            return cls.__as_genomic_feature(q.first())
        return q.first()
    
    @classmethod 
    def select_by_pathogenicity(cls, session, as_genomic_feature=False):
        """returns all strs that are pathogenic"""
        q = session.query(cls).filter(cls.pathogenicity != 0)
        if as_genomic_feature:
            feats = []
            # For each feature...
            for feat in q.all():
                feats.append(
                    cls.__as_genomic_feature(feat)
                )
            return feats
        return q.all() 

    @classmethod
    def select_by_motif(cls, session, motif, compute_rotations=False, as_genomic_feature=False):
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
        if as_genomic_feature:
            feats = []
            # For each feature...
            for feat in q.all():
                feats.append(
                    cls.__as_genomic_feature(feat)
                )
            return feats
        return q.all()  

    @classmethod
    def select_by_uid(cls, session, uid, as_genomic_feature=False):
        q = session.query(cls, Region).\
        join().\
        filter(Region.uid == cls.regionID).\
        filter(cls.uid == uid)
        if as_genomic_feature:
            return cls.__as_genomic_feature(q.first())
        return q.first()

    @classmethod
    def is_unique(cls, session, regionID, sourceID, pathogenicity):
        q = session.query(cls).filter(cls.regionID == regionID, cls.sourceID == sourceID, cls.pathogenicity == pathogenicity)
        q = q.all()
        return len(q) == 0

    @classmethod
    def __as_genomic_feature(self, feat):
        # Define qualifiers
        qualifiers = {
            "uid": feat.ShortTandemRepeat.uid, 
            "regionID": feat.ShortTandemRepeat.regionID, 
            "sourceID": feat.ShortTandemRepeat.sourceID, 
            "motif": feat.ShortTandemRepeat.motif, 
            "pathogenicity": feat.ShortTandemRepeat.pathogenicity         
        }
        return GenomicFeature(
            feat.Region.chrom,
            int(feat.Region.start),
            int(feat.Region.end),
            strand = feat.Region.strand,
            feat_type = "ShortTandemRepeat",
            qualifiers = qualifiers)

    def __repr__(self):
        return "<ShortTandemRepeat(uid={}, regionID={}, sourceID={}, motif={}, pathogencity={})>".format(
            self.uid, self.regionID, self.sourceID, self.motif, self.pathogenicity)
