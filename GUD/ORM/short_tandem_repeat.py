from binning import (
    containing_bins,
    contained_bins,
    assign_bin
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
from .genomicFeatureMixin1 import GFMixin1
from sqlalchemy.ext.declarative import declared_attr


class ShortTandemRepeat(GFMixin1, Base):

    __tablename__ = "short_tandem_repeats"

    # inherits uid, regionID, sourceID
    motif = Column("motif", String(30), nullable=False)
    pathogenicity = Column("pathogenicity", mysql.INTEGER(
        unsigned=True), nullable=False)

    @declared_attr
    def __table_args__(cls):
        return(
            UniqueConstraint(cls.region_id, cls.source_id, cls.pathogenicity),

            Index("ix_str", cls.region_id),
            Index("ix_str_pathogenic", cls.pathogenicity),
            Index("ix_str_motif", cls.motif),

            {
                "mysql_engine": "MyISAM",
                "mysql_charset": "utf8"
            }
        )

    @classmethod
    def select_by_pathogenicity(cls, session, limit, offset):
        """returns all strs that are pathogenic"""
        q = session.query(cls, Region, Source).\
            join().\
            filter(cls.pathogenicity != 0).\
            filter(Region.uid == cls.region_id, Source.uid == cls.source_id)

        return (q.count(), q.offset(offset).limit(limit))

    @classmethod
    def select_by_motif(cls, session, motif, limit, offset, compute_rotations=False):
        """returns all occurences of a certain motif computing 
        rotations if requested"""
        if compute_rotations is True:  # compute the rotations
            motif_len = len(motif)
            motif_temp = motif
            motifs = []
            for i in range(motif_len):
                motif_temp = motif_temp[1:motif_len] + motif_temp[0]
                motifs.append(motif_temp)
        else:
            motifs = [motif]

        q = session.query(cls, Region, Source).\
            join().\
            filter(cls.motif.in_(motifs)).\
            filter(Region.uid == cls.region_id, Source.uid == cls.source_id)

        return (q.count(), q.offset(offset).limit(limit))

    # not in REST API
    @classmethod
    def is_unique(cls, session, regionID, sourceID, pathogenicity):
        q = session.query(cls).filter(cls.region_id == regionID,
                                      cls.source_id == sourceID, cls.pathogenicity == pathogenicity)
        q = q.all()
        return len(q) == 0

    @classmethod
    def as_genomic_feature(self, feat):
        # Define qualifiers
        qualifiers = {
            "uid": feat.ShortTandemRepeat.uid,
            "regionID": feat.ShortTandemRepeat.region_id,
            "sourceID": feat.ShortTandemRepeat.source_id,
            "motif": feat.ShortTandemRepeat.motif,
            "pathogenicity": feat.ShortTandemRepeat.pathogenicity
        }
        return GenomicFeature(
            feat.Region.chrom,
            int(feat.Region.start) + 1,
            int(feat.Region.end),
            strand=feat.Region.strand,
            feat_type="ShortTandemRepeat",
            feat_id="%s_%s" % (self.__tablename__, feat.ShortTandemRepeat.uid),
            qualifiers=qualifiers)
