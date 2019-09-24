from sqlalchemy import (
    Column,
    ForeignKey,
    Integer,
    PrimaryKeyConstraint
)

from sqlalchemy.dialects import mysql
from .base import Base
from .genomic_feature import GenomicFeature
from .region import Region
from .source import Source

from sqlalchemy.ext.declarative import declared_attr


class GFMixin1(object):
    @declared_attr
    def __tablename__(cls):
        return cls.__name__.lower()

    uid = Column("uid", mysql.INTEGER(unsigned=True), primary_key=True)

    @declared_attr
    def region_id(cls):
        return Column("regionID", mysql.INTEGER(unsigned=True), ForeignKey("regions.uid"),
                      nullable=False)

    @declared_attr
    def source_id(cls):
        return Column("sourceID", mysql.INTEGER(unsigned=True), ForeignKey("sources.uid"),
                      nullable=False)

    @classmethod
    def make_query(cls, session):
        q = session.query(cls, Region, Source)\
            .join()\
            .filter(Region.uid == cls.region_id, Source.uid == cls.source_id,)
        return q 

    @classmethod
    def select_by_overlapping_location(cls, session, chrom, start, end):
        """
        Query objects by genomic location.
        """
        bins = Region._compute_bins(start, end)

        q = cls.make_query(session)
        q = q.filter(Region.chrom == chrom,
                    Region.start < end,
                    Region.end > start)\
            .filter(Region.bin.in_(bins))

        return q

    @classmethod
    def select_by_within_location(cls, session, chrom, start, end):
        """
        Query objects by genomic location.
        """
        bins = Region._compute_bins(start, end)

        q = cls.make_query(session)
        q = q.filter(Region.chrom == chrom,
                    Region.start > start,
                    Region.end < end)\
            .filter(Region.bin.in_(bins))

        return q

    @classmethod
    def select_by_exact_location(cls, session, chrom, start, end):
        """
        Query objects by genomic location.
        """
        bins = Region._compute_bins(start, end)

        q = cls.make_query(session)
        q = q.filter(Region.bin.in_(bins))\
            .filter(Region.chrom == chrom,
                    Region.start == start,
                    Region.end == end)

        return q

    @classmethod
    def select_by_location(cls, session, chrom, start, end, location):
        """
        Query objects by genomic location.
        """

        if location == 'exact':
            q = cls.select_by_exact_location(session, chrom, start, end)
        elif location == 'within':
            q = cls.select_by_within_location(session, chrom, start, end)
        elif location == 'overlapping':
            q = cls.select_by_overlapping_location(session, chrom, start, end)
        else:
            return False

        return q

    @classmethod
    def select_by_uids(cls, session, query, uids):
        """
        filter query by uids.
        """
        if (query is None):
            q = cls.make_query(session)
        else: 
            q = query 
        q = q.filter(cls.uid.in_(uids))
        return q

    @classmethod
    def select_by_sources(cls, session, query, sources):
        """
       filter query by sources.
        """
        if(query is None):
            q = cls.make_query(session)
        q = query.filter(Source.name.in_(sources))

        return q

    # for insertion only not in REST
    # implement in child class
    @classmethod
    def is_unique(cls, session):
        return 0

    # for insertion only not in REST
    # implement in child class
    @classmethod
    def select_unique(cls, session):
        return 0

    @classmethod
    def as_genomic_feature(self, feat):
        return GenomicFeature(
            feat.Region.chrom,
            int(feat.Region.start),
            int(feat.Region.end),
            strand=feat.Region.strand,
            feat_type=self.__tablename__,
            feat_id="%s_%s" % (self.__tablename__, self.uid),
            qualifiers=None
        )
