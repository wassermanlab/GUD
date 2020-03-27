from sqlalchemy import (
    Column,
    DateTime,
    Index,
    PrimaryKeyConstraint,
    String,
    UniqueConstraint
)
from sqlalchemy.dialects import mysql
from .base import Base
from datetime import datetime


class Source(Base):

    __tablename__ = "sources"

    uid = Column("uid", mysql.INTEGER(unsigned=True), nullable=False)
    name = Column("name", String(250), nullable=False)
    source_metadata = Column("source_metadata", String(250))
    metadata_descriptor = Column("metadata_descriptor", String(250))
    url = Column("url", String(250))
    insert_date = Column("insert_date", DateTime,
                         default=datetime.utcnow, nullable=False)

    __table_args__ = (
        PrimaryKeyConstraint(uid),
        UniqueConstraint(name, source_metadata, metadata_descriptor, url),
        Index("ix_uid", uid),
        Index("ix_name", name),
        {"mysql_engine": "InnoDB", "mysql_charset": "utf8"}
    )

    @classmethod
    def is_unique(cls, session, name, source_metadata, metadata_descriptor, url):

        q = session.query(cls)\
            .filter(cls.name == name, cls.source_metadata == source_metadata,
                    cls.metadata_descriptor == metadata_descriptor, cls.url == url)

        return len(q.all()) == 0

    @classmethod
    def select_unique(cls, session, name, source_metadata, metadata_descriptor, url):

        q = session.query(cls)\
            .filter(cls.name == name, cls.source_metadata == source_metadata,
                    cls.metadata_descriptor == metadata_descriptor, cls.url == url)

        return q.first()

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

        return "<%s(%s, %s)>" % \
            (
                self.__tablename__,
                "uid={}".format(self.uid),
                "name={}".format(self.name)
            )

    def serialize(self):
        return {
            'uid': self.uid,
            'name': self.name,
            'source_metadata': self.source_metadata, 
            'metadata_descriptor': self.metadata_descriptor, 
            'url': self.url,
            'insert_date': self.insert_date
        }