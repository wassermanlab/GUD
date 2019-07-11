from sqlalchemy import (
    Boolean,
    Column,
    Index,
    PrimaryKeyConstraint,
    String,
    UniqueConstraint
)
from sqlalchemy.dialects import mysql
from sqlalchemy_fulltext import (
    FullText,
    FullTextSearch
)
import sqlalchemy_fulltext.modes \
    as FullTextMode

from .base import Base

class Sample(Base):

    __tablename__ = "samples"

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

    treatment = Column(
        "treatment",
        Boolean,
        nullable=False
    )

    cell_line = Column(
        "cell_line",
        Boolean,
        nullable=False
    )

    cancer = Column(
        "cancer",
        Boolean,
        nullable=False
    )

    __table_args__ = (
        PrimaryKeyConstraint(uid),
        UniqueConstraint(
            name,
            treatment,
            cell_line,
            cancer
        ),
        Index("ix_name", name),
        Index("ix_uid", uid),
        Index(
            "ix_name_fulltext",
            name,
            mysql_prefix="FULLTEXT"
        ),
        Index(
            "ix_treatment_cell_line_cancer",
            treatment,
            cell_line,
            cancer
        ),
        {
            "mysql_engine": "MyISAM",
            "mysql_charset": "utf8"
        }
    )

    @classmethod
    def is_unique(cls, session, name, treatment,
        cell_line, cancer):

        q = session.query(cls)\
            .filter(
                cls.name == name,
                cls.treatment == int(treatment),
                cls.cell_line == int(cell_line),
                cls.cancer == int(cancer)
            )

        return len(q.all()) == 0

    @classmethod
    def select_unique(cls, session, name, treatment,
        cell_line, cancer):

        q = session.query(cls)\
            .filter(
                cls.name == name,
                cls.treatment == int(treatment),
                cls.cell_line == int(cell_line),
                cls.cancer == int(cancer)
            )

        return q.first()

    @classmethod
    def select_by_uid(cls, session, uid):
        """
        Query objects by uid.
        """

        q = session.query(cls).filter(
            cls.uid == uid
        )

        return q.first()

    @classmethod
    def select_by_uids(cls, session, uids=[]):
        """
        Query objects by multiple uids.
        If no uids are provided, return all
        objects.
        """

        q = session.query(cls)

        if uids:
            q = q.filter(cls.uid.in_(uids))

        return q.all()

    @classmethod
    def select_by_name(cls, session, name,
        treatment=False, cell_line=False,
        cancer=False):
        """
        Query objects by sample name. 
        """

        q = session.query(cls)\
            .filter(
                cls.name == name,
                cls.treatment <= int(treatment),
                cls.cell_line <= int(cell_line),
                cls.cancer <= int(cancer)
            )

        return q.all()

    @classmethod
    def select_by_names(cls, session, names=[],
        treatment=False, cell_line=False,
        cancer=False):
        """
        Query objects by multiple sample names.
        If no names are provided, return all
        objects.
        """

        q = session.query(cls)\
            .filter(
                cls.treatment <= int(treatment),
                cls.cell_line <= int(cell_line),
                cls.cancer <= int(cancer)
            )

        if names:
            q = q.filter(cls.name.in_(names))

        return q.all()

    @classmethod
    def select_by_fulltext(
        cls, session, fulltext):
        """
        Query objects by fulltext. 
        """

        class SampleName(FullText, cls):        
            __fulltext_columns__ = list(["name"])

        q = session.query(cls)\
            .filter(
                FullTextSearch(
                    fulltext,
                    SampleName,
                    FullTextMode.NATURAL
                )
            )

        return q.all()

    @classmethod
    def select_by_exp_conditions(cls, session,
        treatment=False, cell_line=False,
        cancer=False):
        """
        Query objects by experimental conditions. 
        """

        q = session.query(cls)\
            .filter(
                cls.treatment <= int(treatment),
                cls.cell_line <= int(cell_line),
                cls.cancer <= int(cancer)
            )

        return q.all()

    def __repr__(self):

        return "<Sample(%s, %s, %s, %s, %s)>" % \
            (
                "uid={}".format(self.uid),
                "name={}".format(self.name),
                "treatment={}".format(self.treatment),
                "cell_line={}".format(self.cell_line),
                "cancer={}".format(self.cancer)
            )