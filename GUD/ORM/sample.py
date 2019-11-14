from fuzzywuzzy import fuzz
from numpy import sqrt
from sqlalchemy import (
    Boolean,
    Column,
    Index,
    Integer,
    PrimaryKeyConstraint,
    String,
    UniqueConstraint
)
from sqlalchemy.dialects import mysql
from sqlalchemy_fulltext import (
    FullText,
    FullTextSearch
)
# import sqlalchemy_fulltext.modes \
#     as FullTextMode

from .base import Base

class Sample(Base):
    # table declerations
    __tablename__ = "samples"

    uid = Column("uid",mysql.INTEGER(unsigned=True),nullable=False)
    name = Column("name",String(250),nullable=False)
    X = Column("X_chroms", mysql.SMALLINT(unsigned=True))
    Y = Column("Y_chroms", mysql.SMALLINT(unsigned=True))
    treatment = Column("treatment", Boolean, nullable=False)
    cell_line = Column("cell_line",Boolean,nullable=False)
    cancer = Column("cancer",Boolean,nullable=False)

    __table_args__ = (
        PrimaryKeyConstraint(uid),
        UniqueConstraint(name, X, Y, treatment, cell_line, cancer),
        Index("ix_name", name),
        Index("ix_uid", uid),
        Index("ix_sex", X, Y),
        Index("ix_treatment_cell_line_cancer", treatment, cell_line, cancer),
        {"mysql_engine": "InnoDB", "mysql_charset": "utf8"}
    )

    @classmethod
    def select_all_samples(cls, session):
        """
        Select all samples
        """

        q = session.query(cls)

        return q

    @classmethod
    def is_unique(cls, session, name, X, Y, treatment,
        cell_line, cancer):

        if X:
            if X.isdigit():
                X = int(X)
        if Y:
            if Y.isdigit():
                Y = int(Y)
    
        q = session.query(cls)\
            .filter(
                cls.name == name,
                cls.X == X,
                cls.Y == Y,
                cls.treatment == int(treatment),
                cls.cell_line == int(cell_line),
                cls.cancer == int(cancer)
            )

        return len(q.all()) == 0

    @classmethod
    def select_unique(cls, session, name, X, Y, treatment,
        cell_line, cancer):

        if X:
            if X.isdigit():
                X = int(X)
        if Y:
            if Y.isdigit():
                Y = int(Y)

        q = session.query(cls)\
            .filter(
                cls.name == name,
                cls.X == X,
                cls.Y == Y,
                cls.treatment == int(treatment),
                cls.cell_line == int(cell_line),
                cls.cancer == int(cancer)
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
    def select_by_fuzzy_string_matching(cls, session, string, threshold=60):
        """
        Returns {Sample}s, and their relevance to the given string, over a
        given threshold.
        """

        samples = []

        for sample in cls.select_all_samples(session).all():
            ratio = fuzz.ratio(string, sample.name)
            partial_ratio = fuzz.partial_ratio(string, sample.name)
            samples.append([sample, sqrt(ratio * partial_ratio)])

        samples.sort(key=lambda x: x[-1], reverse=True)

        return([samples[i] for i in range(len(samples)) if samples[i][1] > threshold])

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

        return "<%s(%s, %s, %s, %s, %s)>" % \
            (
                self.__tablename__,
                "uid={}".format(self.uid),
                "name={}".format(self.name),
                "treatment={}".format(self.treatment),
                "cell_line={}".format(self.cell_line),
                "cancer={}".format(self.cancer)
            )

    def serialize(self):
        return {
            'uid': self.uid,
            'name': self.name,
            'treatment': self.treatment,
            'Y': self.Y,
            'X': self.X,
            'cell_line': self.cell_line,
            'cancer': self.cancer,
            }
