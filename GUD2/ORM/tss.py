from sqlalchemy import (
    Column, Index, PrimaryKeyConstraint, String, ForeignKey,
    UniqueConstraint, Integer, Float, and_, or_
)
from sqlalchemy.dialects import mysql
from GUD2.ORM.region import Region
from GUD2.ORM.source import Source
from GUD2.ORM.sample import Sample
from GUD2.ORM.experiment import Experiment
from GUD2.ORM.gene import Gene 
from GUD2.ORM.base import Base
from binning import containing_bins, contained_bins

class TSS(Base):

    __tablename__ = "tss"

    uid = Column("uid", mysql.INTEGER(unsigned=True))
    regionID = Column("regionID", Integer, ForeignKey("regions.uid"), nullable=False)
    sourceID = Column("sourceID", Integer, ForeignKey("sources.uid"), nullable=False)
    sampleID = Column("sampleID", Integer, ForeignKey("samples.uid"), nullable=False)
    strand = Column("strand", mysql.CHAR(1), nullable=False)
    experimentID = Column("experimentID", Integer, ForeignKey("experiments.uid"), nullable=False)
    gene = Column("gene", String(75), ForeignKey("genes.name2"))
    tss = Column("tss", mysql.INTEGER(unsigned=True))
    avg_tpm = Column("avg_tpm", Float, nullable=False)
    avg_tpm = Column("rel_tpm", Float, nullable=False)

    __table_args__ = (
        PrimaryKeyConstraint(uid),
        UniqueConstraint(regionID, sourceID, sampleID, experimentID),

        Index("ix_tss", regionID), ## query by bin range 
        Index("ix_tss_gene", gene, tss),

        {
            "mysql_engine": "MyISAM",
            "mysql_charset": "utf8"
        }
    )

    @classmethod
    def is_unique(cls, session, regionID, sourceID, sampleID, experimentID):

        q = session.query(cls).filter(
            cls.regionID == regionID,
            cls.sourceID == sourceID,
            cls.sampleID == sampleID,
            cls.experimentID == experimentID
        )

        return len(q.all()) == 0

    @classmethod
    def select_by_sample(cls, session, sample, min_tpm=10.0):

        q = session.query(cls).filter(
            cls.sampleID == sample,
            cls.avg_tpm >= min_tpm
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