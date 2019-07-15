from sqlalchemy import (
    Column,
    Index,
    PrimaryKeyConstraint,
    String,
    UniqueConstraint
)
from sqlalchemy.dialects import mysql

from .base import Base

class Source(Base):

    __tablename__ = "sources"

    uid = Column("uid",mysql.INTEGER(unsigned=True),nullable=False)

    name = Column("name",String(250),nullable=False)

    __table_args__ = (
        PrimaryKeyConstraint(uid),
        UniqueConstraint(name),
        Index("ix_uid", uid),
        Index("ix_name", name),
        {
            "mysql_engine": "MyISAM",
            "mysql_charset": "utf8"
        }
    )

    @classmethod
    def is_unique(cls, session, name):

        q = session.query(cls)\
            .filter(cls.name == name)

        return len(q.all()) == 0

    @classmethod 
    def select_unique(cls, session, name):

        return cls.select_by_name(session, name)

    @classmethod
    def select_by_names(cls, session, names=[]):
        """
        Query objects by multiple source names.
        If no names are provided, return all
        objects.
        """

        q = session.query(cls)

        if names:
            q = q.filter(cls.name.in_(names))

        return q.all()

    @classmethod
    def select_all_sources(cls, session):
        """
        Select all sources
        """

        q = session.query(cls)

        return q

    def __repr__(self):

        return "<Source(%s, %s)>" % \
            (
                "uid={}".format(self.uid),
                "name={}".format(self.name)
            )

    def serialize(self):
        return {
            'uid': self.uid,
            'name': self.name,
            }