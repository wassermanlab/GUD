from sqlalchemy import (
    Column, Index, PrimaryKeyConstraint, String, UniqueConstraint)
from sqlalchemy.dialects import mysql
from .base import Base

class Experiment(Base):
    # table declerations
    __tablename__ = "experiments"
    uid = Column("uid", mysql.INTEGER(unsigned=True), nullable=False)
    name = Column("name", String(250), nullable=False)
    experiment_metadata = Column("experiment_metadata", String(250))
    metadata_descriptor = Column("metadata_descriptor", String(250))
    
    __table_args__ = (
        PrimaryKeyConstraint(uid),
        UniqueConstraint(name, experiment_metadata, metadata_descriptor),
        Index("ix_uid", uid),
        Index("ix_name", name),
        {"mysql_engine": "InnoDB", "mysql_charset": "utf8"}
    )

    # class methods
    @classmethod
    def select_all_experiments(cls, session):
        """
        Select all samples
        """
        q = session.query(cls)
        return q

    @classmethod
    def is_unique(cls, session, name, experiment_metadata, metadata_descriptor):
        """
        returns true if it is unique
        """
        q = session.query(cls).filter(cls.name == name, cls.experiment_metadata == experiment_metadata,
                                      cls.metadata_descriptor == metadata_descriptor)

        return len(q.all()) == 0

    @classmethod
    def select_unique(cls, session, name, experiment_metadata, metadata_descriptor):
        """
        returns selects first by name 
        """
        q = session.query(cls).filter(cls.name == name, cls.experiment_metadata == experiment_metadata,
                                      cls.metadata_descriptor == metadata_descriptor)

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
            ("uid={}".format(self.uid), "name={}".format(self.name))

    def serialize(self):
        return {'uid': self.uid, 'name': self.name, }
