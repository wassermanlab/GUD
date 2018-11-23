from sqlalchemy import (
    Column, Index, PrimaryKeyConstraint, String,
    UniqueConstraint
)
from sqlalchemy.dialects import mysql
from GUD2.ORM.base import Base

class Source(Base):

    __tablename__ = "source"

    uid = Column('uid', mysql.INTEGER(unsigned=True), nullable=False)
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
    def select_by_bin_range(cls, session, name):
        """
        Query objects by name of sample type. 
        """
        q = session.query(cls).filter(cls.name == name, )

        return q.all()

    def __str__(self):
        return "{}".format(self.name)

    def __repr__(self):
        return "<Source(uid={}, name={})>".format(
            self.uid, self.name)
