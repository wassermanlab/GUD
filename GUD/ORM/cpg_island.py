from sqlalchemy import (
    Column,
    Index,
)

from .base import Base
from .region import Region
from .source import Source
from .genomic_feature import GenomicFeature
from .genomicFeatureMixin1 import GFMixin1
from sqlalchemy.ext.declarative import declared_attr

class CpGIsland(GFMixin1, Base):

    __tablename__ = "cpg_islands"
   
    @declared_attr
    def __table_args__(cls):
        return (
        Index("ix_join", cls.source_id, cls.region_id),
        {"mysql_engine": "InnoDB", "mysql_charset": "utf8"}
    )

    #not included in REST API
    @classmethod
    def is_unique(cls, session, regionID, sourceID):
        q = session.query(cls).filter(cls.region_id == regionID,
                                      cls.source_id == sourceID)
        q = q.all()
        return len(q) == 0

    @classmethod
    def as_genomic_feature(self, feat):

        qualifiers = {   
            "uid": feat.CpGIsland.uid,
            "source": feat.Source.name
        }
        return GenomicFeature(
            feat.Region.chrom,
            int(feat.Region.start) ,
            int(feat.Region.end),
            strand = feat.Region.strand,
            feat_type = "CpGIsland",
            feat_id = "%s_%s"%(self.__tablename__, feat.CpGIsland.uid),
            qualifiers = qualifiers)