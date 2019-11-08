from sqlalchemy import (Column, Index, Integer)
from sqlalchemy.dialects import mysql
from .base import Base
from .source import Source
from .genomicFeatureMixin1 import GFMixin1
from sqlalchemy.ext.declarative import declared_attr


class CNV(GFMixin1, Base):
    # table declerations 
    __tablename__ = "copy_number_variants"

    copy_number_change = Column("copy_number_change", Integer, nullable=False)
    clinical_assertion = Column(
        "clinical_assertion", mysql.LONGBLOB, nullable=False)
    clinvar_accession = Column(
        "clinvar_accession", mysql.LONGBLOB, nullable=False)
    dbVar_accession = Column("dbVar_accession", mysql.LONGBLOB, nullable=False)

    @declared_attr
    def __table_args__(cls):
        return (
            Index("ix_join", cls.source_id, cls.region_id),
            {"mysql_engine": "InnoDB", "mysql_charset": "utf8"})

    # class methods
    @classmethod
    def select_by_clinical_assertion(cls, session, query, clinical_assertion):
        """
        add filter to query by clinical assertion 
        """
        q = cls.make_query(session, query)
        q = query.filter(cls.clinical_assertion ==
                         clinical_assertion.encode(encoding='UTF-8'))
        return q

    @classmethod
    def select_by_clinvar_accession(cls, session, query, clinvar_accession):
        """
        filter query by clinvar accession
        """
        q = cls.make_query(session, query)
        accession = ("%" + clinvar_accession + "%").encode(encoding='UTF-8')
        q = query.filter(cls.clinvar_accession.like(accession))
        return q

    @classmethod
    def select_by_dbvar_accession(cls, session, query, dbVar_accession):
        """
        filter query by dbVar_accession
        """
        q = cls.make_query(session, query)
        accession = ("%" + dbVar_accession + "%").encode(encoding='UTF-8')
        q = query.filter(cls.dbVar_accession.like(accession))
        return q

    # not included in REST API
    @classmethod
    def is_unique(cls, session, regionID, sourceID, copy_number_change):
        """
        checks if unique by region, source, and copy number change
        """
        q = session.query(cls).filter(cls.region_id == regionID,
                                      cls.source_id == sourceID,
                                      cls.copy_number_change == copy_number_change)
        q = q.all()
        return len(q) == 0

    @classmethod
    def as_genomic_feature(self, feat):
        """
        extend parent class by adding qualifiers
        """
        qualifiers = {
            "uid": feat.CNV.uid,
            "source": feat.Source.name,
            "copy_number_change": feat.CNV.copy_number_change,
            "clinical_assertion": feat.CNV.clinical_assertion,
            "clinvar_accession": feat.CNV.clinvar_accession,
            "dbVar_accession": feat.CNV.dbVar_accession
        }
        genomic_feature = super().as_genomic_feature(feat)
        genomic_feature.qualifiers = qualifiers
        return genomic_feature
