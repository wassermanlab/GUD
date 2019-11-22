from datetime import datetime
from sqlalchemy import (Column, ForeignKey, String, PrimaryKeyConstraint,
                        UniqueConstraint, Index, DateTime)
from sqlalchemy.dialects import mysql
from .source import Source
from .base import Base
from sqlalchemy.ext.declarative import declared_attr


class Metadata(Base):
    # table declerations
    __tablename__ = "metadata"

    uid = Column("uid", mysql.INTEGER(unsigned=True), nullable=False)
    accession = Column("accession", String(250), nullable=False)
    source_id = Column("sourceID", ForeignKey("sources.uid"), nullable=False)
    insert_date = Column("insert_date", DateTime, default=datetime.utcnow)

    __table_args__ = (
        PrimaryKeyConstraint(uid),
        UniqueConstraint(accession, source_id),   
        {"mysql_engine": "InnoDB", "mysql_charset": "utf8"}
    )

    @classmethod
    def is_unique(cls, session, accession, sourceID):

        q = session.query(cls)\
            .filter(
                cls.accession == accession,
                cls.source_id == sourceID
        )

        return len(q.all()) == 0

    def __repr__(self):
        return "<%s(%s, %s, %s)>" % \
            (self.__tablename__, "accession={}".format(self.accession),
             "sourceID={}".format(self.source_id),
             "insert_date={}".format(self.insert_date))

    def serialize(self):
        return {'accession': self.accession, "sourceID": self.source_id, 
                'insert_date': self.insert_date, }
