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


class GF(object):
    @declared_attr
    def __tablename__(cls):
        return cls.__name__.lower()

    uid = Column("uid", mysql.INTEGER(unsigned=True), primary_key=True)

    @declared_attr
    def region_id(cls):
        return Column("regionID", Integer, ForeignKey("regions.uid"),
                      nullable=False)

    @declared_attr
    def source_id(cls):
        return Column("sourceID", Integer, ForeignKey("sources.uid"),
                      nullable=False)

    @classmethod
    def select_by_location(cls, session, chrom, start, end, limit, offset):
        """
        Query objects by genomic location.
        """
        bins = Region._compute_bins(start, end)

        q = session.query(cls, Region, Source)\
            .join()\
            .filter(Region.uid == cls.region_id, Source.uid == cls.source_id,)\
            .filter(Region.chrom == chrom, 
                    Region.start < end, 
                    Region.end > start)\
            .filter(Region.bin.in_(bins))
            
        return (q.count(), q.offset(offset).limit(limit))

    @classmethod
    def select_by_exact_location(cls, session, chrom, start, end, limit, offset):
        """
        Query objects by genomic location.
        """
        bins = Region._compute_bins(start, end)

        q = session.query(cls, Region, Source)\
            .join()\
            .filter(Region.bin.in_(bins))\
            .filter(Region.uid == cls.region_id, Source.uid == cls.source_id,)\
            .filter(Region.chrom == chrom, 
                    Region.start == start, 
                    Region.end == end)
            
        return (q.count(), q.offset(offset).limit(limit))

    @classmethod
    def select_by_uids(cls, session, uids, limit, offset):
        """
        Query objects by uids.
        """
        q = session.query(cls, Region, Source).\
            join()\
            .filter(Region.uid == cls.region_id, Source.uid == cls.source_id,)\
            .filter(cls.uid.in_(uids))
            
        return (q.count(), q.offset(offset).limit(limit))

    @classmethod
    def select_by_sources(cls, session, sources, limit, offset):
        """
        Query objects by sources.
        """
        q = session.query(cls, Region, Source).\
            join()\
            .filter(Region.uid == cls.region_id, Source.uid == cls.source_id,)\
            .filter(Source.name.in_(sources))
            
        return (q.count(), q.offset(offset).limit(limit))

    @classmethod
    def as_genomic_feature(self, feat):
        return GenomicFeature(
            feat.Region.chrom,
            int(feat.Region.start) + 1,
            int(feat.Region.end),
            strand=feat.Region.strand,
            feat_type=self.__tablename__,
            feat_id="%s_%s" % (self.__tablename__, self.uid),
            qualifiers=None
        )