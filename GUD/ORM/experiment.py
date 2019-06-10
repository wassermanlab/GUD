from sqlalchemy import (
    Column,
    Index,
    PrimaryKeyConstraint,
    String,
    UniqueConstraint
)
from sqlalchemy.dialects import mysql

from .base import Base

class Experiment(Base):

    __tablename__ = "experiments"

    uid = Column(
        "uid",
        mysql.INTEGER(unsigned=True),
        nullable=False
    )

    name = Column(
        "name",
        String(250),
        nullable=False
    )

    __table_args__ = (
        PrimaryKeyConstraint(uid),
        UniqueConstraint(name),
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

        q = session.query(cls)\
            .filter(cls.name == name)

        return q.first()

    @classmethod
    def select_by_name(cls, session, name):
        """
        Query objects by experiment name.
        """

        q = session.query(cls)\
            .filter(cls.name == name)

        return q.first()

    @classmethod
    def select_by_names(cls, session, names=[]):
        """
        Query objects by multiple experiment names.
        If no names are provided, return all
        objects.
        """

        q = session.query(cls)
        
        if names:
            q = q.filter(cls.name.in_(names))

        return q.all()

    def __repr__(self):

        return "<Experiment(%s, %s)>" % \
            (
                "uid={}".format(self.uid),
                "name={}".format(self.name)
            )