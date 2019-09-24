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

    uid = Column("uid", mysql.INTEGER(unsigned=True))

    tssID = Column("tssID", Integer, ForeignKey(
        "transcription_start_sites.uid"), nullable=False)

    sampleID = Column("sampleID", Integer, ForeignKey(
        "samples.uid"), nullable=False)

    avg_expression_level = Column(
        "avg_expression_level", Float, nullable=False)

    __table_args__ = (
        PrimaryKeyConstraint(uid),
        UniqueConstraint(tssID, sampleID),
        Index("ix_tssID", tssID),
        Index("ix_sampleID", sampleID),
        {
            "mysql_engine": "InnoDB",
            "mysql_charset": "utf8"
        }
    )

    def serialize(self, feat):
        return {
            'uid': feat.Expression.uid,
            'tssID': feat.Expression.tssID,
            'sample': feat.Sample.name,
            'avg_tpm': feat.Expression.avg_expression_level,
            }

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
    def select_all(cls, session, query=None):
        """
        """
        if query is None:
            q = session.query(cls, Sample)\
                .join()

        return q

    @classmethod
    def select_by_samples(cls, session, samples, query = None):
        """
        Query objects with a min. expression in 
        one or more sample names.
        If no samples are provided, return all
        objects with a min. expression in any
        sample.
        """
        if query is not None:
            q = query.filter(Sample.name.in_(samples))
        else:
            q = session.query(cls, Sample)\
                .join()\
                .filter(
                    Sample.uid == cls.sampleID
            )\
                .filter(Sample.name.in_(samples))

        return q

    @classmethod
    def select_by_expression(cls, session, min_tpm, max_tpm, query = None):
        """
        Query objects by expression level
        """
        if query is not None:
            q = query.filter(max_tpm >= cls.avg_expression_level >= min_tpm)
        else:
            q = session.query(cls, Sample)\
                .join()\
                .filter(
                    Sample.uid == cls.sampleID
            )\
                .filter(cls.avg_expression_level >= min_tpm)\
                .filter(max_tpm >= cls.avg_expression_level)

        return q

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
