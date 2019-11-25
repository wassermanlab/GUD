from sqlalchemy import (Column, Index, String, UniqueConstraint)
from sqlalchemy.dialects import mysql

from .base import Base
from .region import Region
from .source import Source
from .genomic_feature import GenomicFeature
from .genomicFeatureMixin1 import GFMixin1
from sqlalchemy.ext.declarative import declared_attr


class RepeatMask(GFMixin1, Base):
    # table declerations 
    __tablename__ = "rmsk"

    score = Column("score", mysql.INTEGER, nullable=False)
    name = Column("name", String(75), nullable=False)
    repeat_class = Column("class", String(75), nullable=False)
    family = Column("family", String(75), nullable=False)
    strand = Column("strand", mysql.CHAR(1))
    
    @declared_attr
    def __table_args__(cls):
        return (
        UniqueConstraint(cls.region_id, cls.source_id, cls.name, cls.strand),
        Index("ix_join", cls.region_id, cls.source_id),
        Index("ix_class", cls.repeat_class),
        Index("ix_family", cls.family),
        {"mysql_engine": "InnoDB", "mysql_charset": "utf8"}
    )

    @classmethod
    def is_unique(cls, session, regionID, sourceID, name, strand):
        """
        Checks uniqueness by region source name and strand
        """
        q = session.query(cls).\
            filter(cls.region_id == regionID, cls.source_id == sourceID, 
                   cls.name == name, cls.strand == strand)

        return len(q.all()) == 0

    # # class methods 
    # @classmethod
    # def as_genomic_feature(self, feat):
    #     """
    #     extend parent class by adding qualifiers
    #     """
    #     qualifiers = {
    #         "uid": feat.RepeatMask.uid,
    #         "source": feat.Source.name,
    #         "swScore": feat.RepeatMask.swScore,
    #         "repName": feat.RepeatMask.repName,
    #         "repClass": feat.RepeatMask.repClass,
    #         "repFamily": feat.RepeatMask.repFamily,
    #         "repClass": feat.RepeatMask.repClass
    #     }
    #     genomic_feature = super().as_genomic_feature(feat)
    #     genomic_feature.qualifiers = qualifiers
    #     return genomic_feature
