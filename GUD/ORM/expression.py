from sqlalchemy import (
    and_,
    or_,
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

        q = session.query(cls)\
            .filter(
                cls.tssID == tssID,
                cls.sampleID == sampleID
            )

        return len(q.all()) == 0

    @classmethod
    def select_unique(cls, session, tssID,
        sampleID):

        q = session.query(cls)\
            .filter(
                cls.tssID == tssID,
                cls.sampleID == sampleID
            )

        return q.first()

    @classmethod
    def select_by_sample(cls,
        session, sample, min_tpm=100.0):
        """
        Query objects with a min. expression in 
        sample name.
        """

        q = session.query(cls, Sample)\
            .join()\
            .filter(
                Sample.uid == cls.sampleID,
                TSS.uid == cls.tssID
            )\
            .filter(Sample.name == sample)\
            .filter(cls.avg_expression_level >= min_tpm)

        return q.all()

    @classmethod
    def select_by_samples(cls, session,
        samples=[], min_tpm=100.0):
        """
        Query objects with a min. expression in 
        one or more sample names.
        If no samples are provided, return all
        objects with a min. expression in any
        sample.
        """

        q = session.query(cls, Sample)\
            .join()\
            .filter(
                Sample.uid == cls.sampleID,
                TSS.uid == cls.tssID
            )\
            .filter(cls.avg_expression_level >= min_tpm)

        if samples:
            q = q.filter(Sample.name.in_(samples))

        return q.all()

    def __repr__(self):

        return "<Expression(%s, %s, %s, %s)>" % \
            (
                "uid={}".format(self.uid),
                "tssID={}".format(self.tssID),
                "sampleID={}".format(self.sampleID),
                "avg_expression_level={}".format(
                    self.avg_expression_level
                )
            )