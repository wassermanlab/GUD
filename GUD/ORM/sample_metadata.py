from sqlalchemy import (Column, ForeignKey, String, PrimaryKeyConstraint,
                        UniqueConstraint, Index)
from sqlalchemy.dialects import mysql
from .source import Source
from .base import Base
from sqlalchemy.ext.declarative import declared_attr


class SampleMetadata(Base):
    # table declerations
    __tablename__ = "sample_metadata"

    uid = Column("uid", mysql.INTEGER(unsigned=True), nullable=False)
    gud_name = Column("gud_name", String(250), nullable=False)
    original_name = Column("original_name", String(250), nullable=False)
    sourceID = Column("sourceID", ForeignKey("sources.uid"), nullable=False)

    __table_args__ = (
        PrimaryKeyConstraint(uid),
        UniqueConstraint(original_name, sourceID),
        Index("ix_gud_name", gud_name),
        {"mysql_engine": "InnoDB", "mysql_charset": "utf8"}
    )

    @classmethod
    def select_all_sample_names(cls, session):
        """
        Select all samples
        """
        q = session.query(cls)

        return q

    @classmethod
    def is_unique(cls, session, original_name, sourceID):

        q = session.query(cls)\
            .filter(
                cls.original_name == original_name,
                cls.sourceID == sourceID,
        )

        return len(q.all()) == 0

    @classmethod
    def is_unique(cls, session, original_name, sourceID):

        q = session.query(cls)\
            .filter(
                cls.original_name == original_name,
                cls.sourceID == sourceID,
        )

        return q.first()

    @classmethod
    def select_by_uids(cls, session, uids=[]):
        """
        Query objects by multiple uids.
        If no uids are provided, return all
        objects.
        """

        q = session.query(cls)

        if uids:
            q = q.filter(cls.uid.in_(uids))

        return q.all()

    def serialize(self):
        return {
            'uid': self.uid,
            'gud_name': self.gud_name,
            'original_name': self.original_name,
            'sourceID': self.sourceID,
        }
