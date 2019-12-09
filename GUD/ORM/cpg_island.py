from sqlalchemy import (
    Column,
    Float,
    Index,
    Integer,
    UniqueConstraint
)

from .base import Base
from .region import Region
from .source import Source
from .genomic_feature import GenomicFeature
from .genomicFeatureMixin1 import GFMixin1
from sqlalchemy.ext.declarative import declared_attr

class CpGIsland(GFMixin1, Base):

    __tablename__ = "cpg_islands"

    cpgs = Column("cpgs", Integer, nullable=False)
    gcs = Column("gcs", Integer, nullable=False)
    percent_cpg = Column("percent_cpg", Float, nullable=False)
    percent_gc = Column("percent_gc", Float, nullable=False)
    obsexp_ratio = Column("obsexp_ratio", Float, nullable=False)

    @declared_attr
    def __table_args__(cls):
        return (
        UniqueConstraint(cls.region_id, cls.source_id),
        Index("ix_join", cls.region_id, cls.source_id),
        {"mysql_engine": "InnoDB", "mysql_charset": "utf8"}
    )

  
    @classmethod
    def as_genomic_feature(self, feat):

        qualifiers = {
            "uid": feat.CpGIsland.uid,
            "cpgs": feat.CpGIsland.cpgs,
            "gcs": feat.CpGIsland.gcs,
            "percent_cpg": feat.CpGIsland.percent_cpg,
            "percent_gc": feat.CpGIsland.percent_gc,
            "obsexp_ratio": feat.CpGIsland.obsexp_ratio,
            "source": feat.sourceName,
        }

        genomic_feature = super().as_genomic_feature(feat)
        genomic_feature.qualifiers = qualifiers
        return genomic_feature