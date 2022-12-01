from sqlalchemy import (Column, Index, String, UniqueConstraint)
from sqlalchemy.dialects import mysql
from .genomicFeatureMixin1 import GFMixin1 # replace with GFMixin2 if necissary 
from .genomic_feature import GenomicFeature
from .base import Base
from .region import Region
from .source import Source
from sqlalchemy.ext.declarative import declared_attr

class FeatureName(GFMixin1, Base):                                              # replace "FeatureName"
    # table declerations
    __tablename__ = "table_name"                                                # replace "table_name"

    # add additional columns                                                 
    # string_column = Column("string_column_name", String(75), nullable=False)
    # int_column = Column("int_column_name", mysql.INTEGER(unsigned=True), nullable=False)

    @declared_attr
    def __table_args__(cls):
        return (
            # add unique constrains if one exists
            # UniqueConstraint(cls.region_id, cls.unique_column, cls.source_id),
            # query by bin range
            Index("ix_source_join", cls.source_id, cls.region_id),
            # add any additional relevent indexes
            {"mysql_engine": "InnoDB", "mysql_charset": "utf8", }
        )

    # add additional class methods class methods
    # @classmethod
    # def select_by_additional_column(cls, session, query, columns=[]):
    #     """
    #     Query objects by multiple column features
    #     """
    #     q = cls.make_query(session, query)  # first make query if one doesn't exists
    #     q = q.filter(cls.name2.in_(columns))  # filter as query is required
    #     return q

    # override this method in parent class if there is a unique condition
    # @classmethod
    # def is_unique(cls, session, regionID, sourceID):
    #     """
    #     Checks uniqueness by region source and name
    #     """
    #     q = session.query(cls).filter(cls.region_id == regionID, cls.source_id == sourceID)
    #     return len(q.all()) == 0

    @classmethod
    def as_genomic_feature(self, feat):
        """
        extend parent class by adding qualifiers
        """
        # Define qualifiers
        qualifiers = {
            # "uid": feat.Gene.uid,             # add all qualifiers to feature (any additional columns of this feature)
        }
        genomic_feature = super().as_genomic_feature(feat)
        genomic_feature.qualifiers = qualifiers
        return genomic_feature
