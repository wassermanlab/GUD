from sqlalchemy import (
    Column, Index, PrimaryKeyConstraint, String,
    UniqueConstraint, CheckConstraint, Boolean
)
from sqlalchemy.dialects import mysql
from sqlalchemy.ext.declarative import declarative_base
Base = declarative_base()


class Sample(Base):

    __tablename__ = "sample"

    uid = Column('uid', mysql.INTEGER(unsigned=True), nullable=False)
    name = Column("name", String(250), nullable=False)
    treatment = Column("treatment", Boolean, nullable=False)
    cell_line = Column("cell_line", Boolean, nullable=False)
    cancer = Column("cancer", Boolean, nullable=False)

    __table_args__ = (
        PrimaryKeyConstraint(uid),
        UniqueConstraint(name, treatment, cell_line, cancer),

        Index("ix_tss", name),

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
        return "{}\t{}\t{}\t{}".format(self.name, self.treatment,
                                       self.cell_line, self.cancer)

    def __repr__(self):
        return "<Sample(uid={}, name={}, treatment={}, cell_line={}, cancer={})>".format(
            self.uid, self.name, self.treatment, self.cell_line, self.cancer)
