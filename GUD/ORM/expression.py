from sqlalchemy import (
    Column,
    Float,
    ForeignKey,
    Index,
    Integer,
    PrimaryKeyConstraint,
    UniqueConstraint
)
from sqlalchemy.dialects import mysql

from .base import Base
from .sample import Sample
from .tss import TSS

class Expression(Base):

    __tablename__ = "expression"

    uid = Column(
        "uid",
        mysql.INTEGER(unsigned=True)
    )

    tssID = Column(
        "tssID",
        Integer,
        ForeignKey("transcription_start_sites.uid"),
        nullable=False
    )

    sampleID = Column(
        "sampleID",
        Integer,
        ForeignKey("samples.uid"),
        nullable=False
    )

    avg_expression_level = Column(
        "avg_expression_level",
        Float,
        nullable=False
    )

    __table_args__ = (
        PrimaryKeyConstraint(uid),
        UniqueConstraint(tssID, sampleID),
        Index("ix_tssID", tssID),
        Index("ix_sampleID", sampleID),
        {
            "mysql_engine": "MyISAM",
            "mysql_charset": "utf8"
        }
    )

    @classmethod
    def is_unique(cls, session, tssID, sampleID):

        q = session.query(cls).\
            filter(
                cls.tssID == tssID,
                cls.sampleID == sampleID
            )

        return len(q.all()) == 0

    @classmethod
    def select_unique(cls, session, tssID, sampleID):

        q = session.query(cls).\
            filter(
                cls.tssID == tssID,
                cls.sampleID == sampleID
            )

        return q.first()

    
#    @classmethod
#    def select_by_sample(cls, session, sample, min_tpm=10.0):
#
#        q = session.query(cls).\
#            filter(
#                cls.sampleID == sample,
#                cls.avg_tpm >= min_tpm
#            )
#
#        return q.first()
#
#    @classmethod
#    def select_by_samples(cls, session, sample=[], min_tpm=10.0):
#
#        q = session.query(cls).\
#            filter(cls.avg_tpm >= min_tpm).\
#            filter(cls.sampleID.in_(sample))
#
#        return q.all()
#
#    @classmethod
#    def select_by_tss(cls, session, tss, sample=[]):
#        """
#        Query objects by TSS id.
#        """
#
#        q = session.query(cls).filter(cls.tssID == tssID)
#
#        if gene and tss:
#            q = q.filter(cls.gene == gene, cls.tss == tss)
#
#        if sample:
#            q = q.filter(cls.sampleID.in_(sample))
#
#        return q.all()
#
#    @classmethod
#    def select_by_multiple_tss(cls, session, tss=[], sample=[]):
#        """
#        Query objects by multiple TSS ids. If no TSS is
#        provided, return all objects.
#        """ 
#
#        q = session.query(cls)
#
#        if tss:
#            # Initialize
#            ands = []
#            # For each gene, TSS pair...
#            for i, j in tss:
#                ands.append(and_(cls.gene == i,
#                    cls.tss == j))
#            q = q.filter(or_(*ands))
#
#        if sample:
#            q = q.filter(cls.sampleID.in_(sample))
#
#        return q.all()
#
#    @classmethod
#    def get_all_samples(cls, session):
#        """
#        Query all TSS objects in the database and return
#        theirsamples (sampleID field).
#        """
#        samples = session.query(cls.sampleID).distinct().all()
#
#        return [s[0] for s in samples]