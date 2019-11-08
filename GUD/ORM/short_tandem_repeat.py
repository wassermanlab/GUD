from sqlalchemy import (Column, Index, String, UniqueConstraint)
from sqlalchemy.dialects import mysql
from .base import Base
from .region import Region
from .source import Source
from .genomicFeatureMixin1 import GFMixin1
from sqlalchemy.ext.declarative import declared_attr


class ShortTandemRepeat(GFMixin1, Base):
    # table declerations
    __tablename__ = "short_tandem_repeats"

    # inherits uid, regionID, sourceID
    motif = Column("motif", String(30), nullable=False)
    pathogenicity = Column("pathogenicity", mysql.INTEGER(
        unsigned=True), nullable=False)

    @declared_attr
    def __table_args__(cls):
        return (
        UniqueConstraint(cls.region_id, cls.source_id, cls.pathogenicity),
        Index("ix_join", cls.region_id, cls.source_id),
        Index("ix_str_pathogenic", cls.pathogenicity),
        Index("ix_str_motif", cls.motif),
        {"mysql_engine": "InnoDB", "mysql_charset": "utf8"}
    )

    #class methods
    @classmethod
    def select_by_pathogenicity(cls, session, query):
        """returns all strs that are pathogenic"""
        q = cls.make_query(session, query)
        q = q.filter(cls.pathogenicity != 0)
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
        q = cls.make_query(session, query)
        q = q.filter(cls.motif.in_(motifs))
        return q

    # not in REST API
    @classmethod
    def is_unique(cls, session, regionID, sourceID, pathogenicity):
        """check uniqueness by region,  source, pathogenicity"""
        q = session.query(cls).filter(cls.region_id == regionID,
                                      cls.source_id == sourceID, 
                                      cls.pathogenicity == pathogenicity)
        q = q.all()
        return len(q) == 0

    @classmethod
    def as_genomic_feature(self, feat):
        """
        extend parent class by adding qualifiers
        """
        qualifiers = {
            "uid": feat.ShortTandemRepeat.uid,
            "source": feat.Source.name,
            "motif": feat.ShortTandemRepeat.motif,
            "pathogenicity": feat.ShortTandemRepeat.pathogenicity
        }
        genomic_feature = super().as_genomic_feature(feat)
        genomic_feature.qualifiers = qualifiers
        return genomic_feature

