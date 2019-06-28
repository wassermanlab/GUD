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

from .base import Base
from .genomic_feature import GenomicFeature
from .region import Region
from .source import Source

class Gene(Base):

    __tablename__ = "genes"

    uid = Column(
        "uid",
        mysql.INTEGER(unsigned=True)
    )

    regionID = Column(
        "regionID",
        Integer,
        ForeignKey("regions.uid"),
        nullable=False
    )

    sourceID = Column(
        "sourceID",
        Integer,
        ForeignKey("sources.uid"),
        nullable=False
    )

    name = Column(
        "name",
        String(75),
        nullable=False
    )

    name2 = Column(
        "name2",
        String(75),
        nullable=False
    )

    cdsStart = Column(
        "cdsStart",
        mysql.INTEGER(unsigned=True),
        nullable=False
    )

    cdsEnd = Column(
        "cdsEnd",
        mysql.INTEGER(unsigned=True),
        nullable=False
    )

    exonStarts = Column(
        "exonStarts",
        mysql.LONGBLOB,
        nullable=False
    )

    exonEnds = Column(
        "exonEnds",
        mysql.LONGBLOB,
        nullable=False
    )

    __table_args__ = (
        PrimaryKeyConstraint(uid),
        UniqueConstraint(
            regionID,
            sourceID,
            name
        ),
        Index("ix_regionID", regionID), # query by bin range
        Index("ix_name", name),
        Index("ix_name2", name2),
        {
            "mysql_engine": "MyISAM",
            "mysql_charset": "utf8"
        }
    )

    @classmethod
    def is_unique(cls, session, regionID, sourceID, name):

        q = session.query(cls)\
            .filter(
                cls.regionID == regionID,
                cls.sourceID == sourceID,
                cls.name == name,
            )

        return len(q.all()) == 0

    @classmethod
    def select_unique(cls, session, regionID, sourceID, name):

        q = session.query(cls)\
            .filter(
                cls.regionID == regionID,
                cls.sourceID == sourceID,
                cls.name == name,
            )

        return q.first()

    @classmethod
    def select_by_location(cls, session, chrom, start, end, as_genomic_feature=False):
        """
        Query objects by genomic location.
        """

        bins = Region._compute_bins(start, end)

        q = session.query(cls, Region, Source)\
            .join()\
            .filter(
                Region.uid == cls.regionID,
                Source.uid == cls.sourceID,
            )\
            .filter(
                Region.chrom == chrom,
                Region.start < end,
                Region.end > start,
            )\
            .filter(Region.bin.in_(bins))

        if as_genomic_feature:

            feats = []

            # For each feature...
            for feat in q.all():
                feats.append(
                    cls.__as_genomic_feature(feat)
                )
            
            return feats
    
        return q.all()

    @classmethod
    def select_by_name(cls, session, name, as_genomic_feature=False):
        """
        Query objects by gene symbol.
        """

        q = session.query(cls, Region, Source)\
            .join()\
            .filter(
                Region.uid == cls.regionID,
                Source.uid == cls.sourceID,
            )\
            .filter(cls.name2 == name)

        if as_genomic_feature:

            feats = []

            # For each feature...
            for feat in q.all():
                feats.append(
                    cls.__as_genomic_feature(feat)
                )

            return feats

        return q.all()

    @classmethod
    def select_by_names(cls, session, names=[], as_genomic_feature=False):
        """
        Query objects by multiple gene symbols.
        If no genes are provided, return all
        objects.
        """

        q = session.query(cls, Region, Source)\
            .join()\
            .filter(
                Region.uid == cls.regionID,
                Source.uid == cls.sourceID,
            )\

        if names:
            q = q.filter(cls.name2.in_(names))

        if as_genomic_feature:

            feats = []

            # For each feature...
            for feat in q.all():
                feats.append(
                    cls.__as_genomic_feature(feat)
                )

            return feats

        return q.all()

    @classmethod
    def select_by_uid(cls, session, uid, as_genomic_feature=False):
        """
        Query objects by uid.
        """

        q = session.query(cls, Region, Source).\
            join().\
            filter(cls.uid == uid,
                Region.uid == cls.regionID,
                Source.uid == cls.sourceID,
            )

        if as_genomic_feature:
            return cls.__as_genomic_feature(
                q.first()
            )

        return q.first()

#
#    @classmethod
#    def select_by_uids(cls, session, uids=[],
#        as_genomic_feature=False):
#        """
#        Query objects by multiple uids.
#        If no uids are provided, return all
#        objects.
#        """
#
#        q = session.query(cls, Region, Source).\
#            join()
#        
#        if uids:
#            q = q.filter(cls.uid.in_(uids))
#
#        if as_genomic_feature:
#
#            feats = []
#
#            # For each feature...
#            for feat in q.all():
#                feats.append(
#                    cls.__as_genomic_feature(feat)
#                )
#
#            return feats
#
#        return q.all()

    @classmethod
    def get_all_gene_symbols(cls, session):
        """
        Return the gene symbol (name2 field) of all
        objects.
        """

        q = session.query(cls.name2).distinct().all()

        return [g[0] for g in q]

    @classmethod
    def __as_genomic_feature(self, feat):

        # Note for Oriol:
        # Move this within your code!!!

#        # Initialize
#        exonStarts = []
#        exonEnds = []
#
#        # For each exon start...
#        for i in str(feat.Gene.exonStarts).split(","):
#            if i.isdigit():
#                exonStarts.append(int(i))
#        
#        # For each exon end...
#        for i in str(feat.Gene.exonEnds).split(","):
#            if i.isdigit():
#                exonStarts.append(int(i))

        # Define qualifiers
        qualifiers = {
            "uid": feat.Gene.uid,
            "regionID": feat.Gene.regionID,
            "sourceID": feat.Gene.sourceID,
            "source": feat.Source.name,
            "name": feat.Gene.name,
            "name2": feat.Gene.name2,
            "cdsStart": int(feat.Gene.cdsStart),
            "cdsEnd": int(feat.Gene.cdsEnd),
            "exonStarts": feat.Gene.exonStarts,
            "exonEnds": feat.Gene.exonEnds,      
        }

        return GenomicFeature(
            feat.Region.chrom,
            int(feat.Region.start),
            int(feat.Region.end),
            strand = feat.Region.strand,
            feat_type = self.__tablename__,
            feat_id = "%s_%s" % \
                (
                    self.__tablename__,
                    feat.Gene.uid
                ),
            qualifiers = qualifiers
        )

    def __repr__(self):

        return "<%s(%s, %s, %s, %s, %s)>" % \
            (   
                self.__tablename__,
                "uid={}".format(self.uid),
                "regionID={}".format(self.regionID),
                "sourceID={}".format(self.sourceID),
                "name={}".format(self.name),
                "name2={}".format(self.name2),
            )