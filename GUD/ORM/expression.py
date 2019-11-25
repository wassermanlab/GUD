from sqlalchemy import (Column, Index, ForeignKey, PrimaryKeyConstraint, UniqueConstraint)
from .base import Base
from .sample import Sample
from .tss import TSS
from sqlalchemy.dialects import mysql


class Expression(Base):

    __tablename__ = "expression"
    uid = Column("uid", mysql.INTEGER(unsigned=True))
    expression_level = Column("expression_level", mysql.FLOAT, nullable=False)
    tss_id = Column("tssID", mysql.INTEGER, ForeignKey("transcription_start_sites.uid"),
                    nullable=False)
    sample_id = Column("sampleID", mysql.INTEGER, ForeignKey("samples.uid"),
                       nullable=False)
 
    __table_args__ = (
        PrimaryKeyConstraint(uid),
        UniqueConstraint(tss_id, sample_id),
        Index("ix_tssID", tss_id),
        Index("ix_sampleID", sample_id),
        {"mysql_engine": "InnoDB", "mysql_charset": "utf8"}
    )

    # class methods
    @classmethod
    def is_unique(cls, session, tssID, sampleID):

        q = session.query(cls).filter(cls.tss_id == tssID,
                                      cls.sample_id == sampleID)

        return len(q.all()) == 0

    @classmethod
    def select_unique(cls, session, tssID, sampleID):

        q = session.query(cls).filter(cls.tss_id == tssID,
                                      cls.sample_id == sampleID)

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
    def select_by_expression(cls, session, min_tpm, max_tpm = None, query = None):
        """
        Query objects by expression level
        """
        if query is not None:
            q = query.filter(cls.expression_level >= min_tpm)
            if max_tpm is not None:
                q = q.filter(max_tpm >= cls.expression_level)
        else:
            q = session.query(cls, Sample)\
                .join()\
                .filter(
                    Sample.uid == cls.sampleID
                )\
                .filter(
                    cls.expression_level >= min_tpm
                )
            if max_tpm is not None:
                q = q.filter(max_tpm >= cls.expression_level)

        return q

    def serialize(self, feat):
        return {
            'uid': feat.Expression.uid,
            'tssID': feat.Expression.tss_id,
            'sample': feat.Sample.name,
            'expression_level': feat.Expression.expression_level,
            }

    def __repr__(self):

        return "<Expression(%s, %s, %s, %s)>" % \
            (
                "uid={}".format(self.uid),
                "tssID={}".format(self.tssID),
                "sampleID={}".format(self.sampleID),
                "expression_level={}".format(
                    self.expression_level
                )
            )
