from sqlalchemy import (
    and_,
    or_,
    Column,
    Index,
    Integer,
    PrimaryKeyConstraint,
    String,
    ForeignKey,
    UniqueConstraint
)
from sqlalchemy.dialects import mysql

from .base import Base
from .experiment import Experiment
from .gene import Gene 
from .region import Region
from .source import Source

class TSS(Base):

    __tablename__ = "transcription_start_sites"

    uid = Column("uid", mysql.INTEGER(unsigned=True))
    regionID = Column("regionID", Integer, ForeignKey("regions.uid"), nullable=False)
    sourceID = Column("sourceID", Integer, ForeignKey("sources.uid"), nullable=False)
    experimentID = Column("experimentID", Integer, ForeignKey("experiments.uid"), nullable=False)
    gene = Column("gene", String(75), ForeignKey("genes.name2"))
    tss = Column("tss", mysql.INTEGER(unsigned=True))
    strand = Column("strand", mysql.CHAR(1), nullable=False)
    samples = Column("sampleIDs", mysql.LONGBLOB, nullable=False)
    samples = Column("avg_tpms", mysql.LONGBLOB, nullable=False)

    __table_args__ = (
        PrimaryKeyConstraint(uid),
        UniqueConstraint(regionID, sourceID, experimentID),

        Index("ix_tss", regionID), # query by bin range 
        Index("ix_tss_gene", gene, tss),

        {
            "mysql_engine": "MyISAM",
            "mysql_charset": "utf8"
        }
    )

    @classmethod
    def is_unique(cls, session, regionID, sourceID, sampleID,
        experimentID):

        q = session.query(cls).filter(
            cls.regionID == regionID,
            cls.sourceID == sourceID,
            cls.sampleID == sampleID,
            cls.experimentID == experimentID
        )

        return len(q.all()) == 0

    @classmethod
    def select_by_sample(cls, session, sample, min_tpm=10.0,
        rel_tpm=0.0):

        q = session.query(cls).filter(
            cls.sampleID == sample,
            cls.avg_tpm >= min_tpm,
            cls.rel_tpm >= rel_tpm
        )

        return q.first()

    @classmethod
    def select_by_samples(cls, session, sample=[], min_tpm=10.0):

        q = session.query(cls).filter(cls.sampleID.in_(sample)).\
            filter(cls.avg_tpm >= min_tpm)

        return q.all()

    @classmethod
    def select_by_tss(cls, session, gene, tss, sample=[]):
        """Query objects by TSS (i.e. gene + tss). If no TSS is
        provided, query all TSSs.
        """
        q = session.query(cls)

        if gene and tss:
            q = q.filter(cls.gene == gene, cls.tss == tss)

        if sample:
            q = q.filter(cls.sampleID.in_(sample))

        return q.all()

    @classmethod
    def select_by_multiple_tss(cls, session, tss=[], sample=[]):
        """Query objects by list of TSSs. If no TSS are provided,
        query all TSSs. Provide TSSs as a {list} of {lists}/{tuples}
        of length 2 in the form gene, tss.
        """    
        q = session.query(cls)

        if tss:
            # Initialize
            ands = []
            # For each gene, TSS pair...
            for i, j in tss:
                ands.append(and_(cls.gene == i,
                    cls.tss == j))
            q = q.filter(or_(*ands))

        if sample:
            q = q.filter(cls.sampleID.in_(sample))

        return q.all()

    @classmethod
    def get_all_samples(cls, session):
        """Query all TSS objects in the database and return their
        samples (sampleID field).
        """
        samples = session.query(cls.sampleID).distinct().all()

        return [s[0] for s in samples]