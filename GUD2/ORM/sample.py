from sqlalchemy import (
    Column, Index, PrimaryKeyConstraint, String,
    UniqueConstraint, CheckConstraint, Boolean
)
from sqlalchemy.dialects import mysql
from GUD2.ORM.base import Base

class Sample(Base):

    __tablename__ = "samples"

    uid = Column("uid", mysql.INTEGER(unsigned=True), nullable=False)
    name = Column("name", String(250), nullable=False)
    treatment = Column("treatment", Boolean, nullable=False)
    cell_line = Column("cell_line", Boolean, nullable=False)
    cancer = Column("cancer", Boolean, nullable=False)

    __table_args__ = (
        PrimaryKeyConstraint(uid),
        UniqueConstraint(name, treatment, cell_line, cancer),

        Index("ix_sample", name),

        {
            "mysql_engine": "MyISAM",
            "mysql_charset": "utf8"
        }
    )

    @classmethod
    def select_by_name(cls, session, name, treatment=False,
        cell_line=False, cancer=False):
        """Query objects by name of sample type. 
        """
        q = session.query(cls).filter(cls.name == name).\
            filter(
                cls.treatment <= int(treatment),
                cls.cell_line <= int(cell_line),
                cls.cancer <= int(cancer)
            )

        return q.first()

    @classmethod
    def select_by_names(cls, session, names, treatment=False,
        cell_line=False, cancer=False):
        """Query objects by names of sample type. 
        """
        q = session.query(cls).filter(cls.name.in_(names)).\
            filter(
                cls.treatment <= int(treatment),
                cls.cell_line <= int(cell_line),
                cls.cancer <= int(cancer)
            )

        return q.all()

    @classmethod
    def select_by_exp_conditions(cls, session, treatment=False,
        cell_line=False, cancer=False):
        """Query objects by experimental conditions. 
        """
        q = session.query(cls).filter(
            cls.treatment <= int(treatment),
            cls.cell_line <= int(cell_line),
            cls.cancer <= int(cancer)
        )

        return q.all()

    @classmethod 
    def select_by_exact_sample(cls, session, name, treatment,
        cell_line, cancer):

        q = session.query(cls).filter(
            cls.name == name,
            cls.treatment == treatment,
            cls.cell_line == cell_line,
            cls.cancer == cancer
        )

        return q.first()

    def __str__(self):
        return "{}\t{}\t{}\t{}".format(self.name, self.treatment,
            self.cell_line, self.cancer)

    def __repr__(self):
        return "<Sample(uid={}, name={}, treatment={}, cell_line={}, cancer={})>".format(
            self.uid, self.name, self.treatment, self.cell_line, self.cancer)