from sqlalchemy import (
    Column,
    ForeignKey,
    Index,
    Integer,
    PrimaryKeyConstraint,
    String,
    UniqueConstraint
)
from sqlalchemy.dialects import mysql

from .gene_feature import GeneFeature
from .genomic_feature import GenomicFeature
from .base import Base
from .region import Region
from .source import Source

from sqlalchemy.ext.declarative import declared_attr


class Gene(GeneFeature, Base):

    __tablename__ = "genes"
    # inherits uid, regionID, sourceID

    name = Column("name", String(75), nullable=False)
    name2 = Column("name2", String(75), nullable=False)
    cdsStart = Column("cdsStart", mysql.INTEGER(unsigned=True), nullable=False)
    cdsEnd = Column("cdsEnd", mysql.INTEGER(unsigned=True), nullable=False)
    exonStarts = Column("exonStarts", mysql.LONGBLOB, nullable=False)
    exonEnds = Column("exonEnds", mysql.LONGBLOB, nullable=False)

    @declared_attr
    def __table_args__(cls):
        return (
            UniqueConstraint(
                cls.region_id,
                cls.name,
                cls.source_id
            ),
            # query by bin range
            Index("ix_regionID", cls.region_id),
            Index("ix_name", cls.name),
            Index("ix_name2", cls.name2),
            {
                "mysql_engine": "MyISAM",
                "mysql_charset": "utf8",
            }
        )

    @classmethod
    def select_by_names(cls, session, names=[]):
        """
        Query objects by multiple gene symbols.
        If no genes are provided, return all
        objects.
        """
        q = session.query(cls, Region, Source)\
            .join()\
            .filter(Region.uid == cls.regionID, 
                    Source.uid == cls.sourceID,)\

        if names:
            q = q.filter(cls.name2.in_(names))

        return q.all()

    @classmethod
    def get_all_gene_symbols(cls, session):
        """
        Return the gene symbol (name2 field) of all objects.
        """
        q = session.query(cls.name2).distinct().all()
        return [g[0] for g in q]

    # for insertion only not in REST
    @classmethod
    def is_unique(cls, session, regionID, name, sourceID):

        q = session.query(cls)\
            .filter(cls.regionID == regionID, 
                    cls.name == name, 
                    cls.sourceID == sourceID)

        return len(q.all()) == 0

    # for insertion only not in REST
    @classmethod
    def select_unique(cls, session, regionID,
                      name, sourceID):

        q = session.query(cls)\
            .filter(cls.regionID == regionID, 
                    cls.name == name, 
                    cls.sourceID == sourceID)

        return q.first()

    @classmethod
    def as_genomic_feature(self, feat):
        # Initialize
        exonStarts = []
        exonEnds = []

        # For each exon start...
        for i in str(feat.Gene.exonStarts).split(","):
            if i.isdigit():
                exonStarts.append(int(i))

        # For each exon end...
        for i in str(feat.Gene.exonEnds).split(","):
            if i.isdigit():
                exonStarts.append(int(i))

        # Define qualifiers
        qualifiers = {
            "uid": feat.Gene.uid,
            "regionID": feat.Gene.region_id,
            "name": feat.Gene.name,
            "name2": feat.Gene.name2,
            "cdsStart": int(feat.Gene.cdsStart),
            "cdsEnd": int(feat.Gene.cdsEnd),
            "exonStarts": feat.Gene.exonStarts,
            "exonEnds": feat.Gene.exonEnds,
            "sourceID": feat.Gene.source_id,
            "source": feat.Source.name,
        }

        return GenomicFeature(
            feat.Region.chrom,
            int(feat.Region.start),
            int(feat.Region.end),
            strand=feat.Region.strand,
            feat_type="Gene",
            feat_id="%s_%s" % (self.__tablename__, feat.Gene.uid),
            qualifiers=qualifiers
        )
