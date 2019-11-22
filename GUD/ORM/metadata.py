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
                cls.sourceID == sourceID
        )

        return len(q.all()) == 0
