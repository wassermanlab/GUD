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

    uid = Column("uid", mysql.INTEGER(unsigned=True), nullable=False)
    name = Column("name", String(250), nullable=False)

    __table_args__ = (
        PrimaryKeyConstraint(uid),
        UniqueConstraint(name),
        Index("ix_source", name),
        {
            "mysql_engine": "MyISAM",
            "mysql_charset": "utf8"
        }
    )

    @classmethod
    def select_by_name(cls, session, name):
        """
        Query objects by source name. 
        """

        q = session.query(cls).filter(cls.name == name)

        return q.first()

    @classmethod 
    def select_unique(cls, session, name):

        return cls.select_by_name(session, name)

    def __str__(self):

        return "{}".format(self.name)

    def __repr__(self):

        return "<Source(%s, %s)>" % \
            (
                "uid={}".format(self.uid),
                "name={}".format(self.name)
            )