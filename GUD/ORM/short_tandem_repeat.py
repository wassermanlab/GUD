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

    __table_args__ = (
        UniqueConstraint(region_id, source_id, pathogenicity),
        Index("ix_join", region_id, source_id),
        Index("ix_str_pathogenic", pathogenicity),
        Index("ix_str_motif", motif),
        {"mysql_engine": "InnoDB", "mysql_charset": "utf8"}
    )

    @classmethod
    def select_by_pathogenicity(cls, session):
        """returns all strs that are pathogenic"""
        q = session.query(cls, Region, Source).\
            join().\
            filter(cls.pathogenicity != 0).\
            filter(Region.uid == cls.region_id, Source.uid == cls.source_id)

        return q

    @classmethod
    def select_by_motif(cls, session, motif, query, compute_rotations=False):
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

        q = query.filter(cls.motif.in_(motifs))

        return q

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
            "source": feat.Source.name,
            "motif": feat.ShortTandemRepeat.motif,
            "pathogenicity": feat.ShortTandemRepeat.pathogenicity
        }
        return GenomicFeature(
            feat.Region.chrom,
            int(feat.Region.start),
            int(feat.Region.end),
            strand=feat.Region.strand,
            feat_type="ShortTandemRepeat",
            feat_id="%s_%s" % (self.__tablename__, feat.ShortTandemRepeat.uid),
            qualifiers=qualifiers)
