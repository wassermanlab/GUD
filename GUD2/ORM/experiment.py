from sqlalchemy import (
    Column,
    Index,
    PrimaryKeyConstraint,
    String,
    UniqueConstraint
)
from sqlalchemy.dialects import mysql

from GUD2.ORM.base import Base

class Experiment(Base):

    __tablename__ = "experiments"

    uid = Column("uid", mysql.INTEGER(unsigned=True), nullable=False)
    name = Column("name", String(250), nullable=False)

    __table_args__ = (
        PrimaryKeyConstraint(uid),
        UniqueConstraint(name),
        Index("ix_experiment", name),
        {
            "mysql_engine": "MyISAM",
            "mysql_charset": "utf8"
        }
    )

    @classmethod
    def select_by_name(cls, session, name):
        """
        Query objects by experiment name.
        """

        q = session.query(cls).filter(cls.name == name)

        return q.first()

    @classmethod
    def select_by_names(cls, session, names=[]):
        """
        Query objects by multiple experiment names. If no
        experiment names are provided, return all objects.
        """

        q = session.query(cls).filter(cls.name.in_(names))

        return q.all()

    def __str__(self):
        return "{}".format(self.name)

    def __repr__(self):
        return "<Experiment(%s, %s)>" % \
            (
                "uid={}".format(self.uid),
                "name={}".format(self.name)
            )