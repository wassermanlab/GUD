from sqlalchemy import (Column, ForeignKey)
from sqlalchemy.dialects import mysql
from .genomic_feature import GenomicFeature
from .region import Region
from .source import Source
from sqlalchemy.ext.declarative import declared_attr


class GFMixin1(object):
    # table declaration
    @declared_attr
    def __tablename__(cls):
        return cls.__name__.lower()

    uid = Column("uid", mysql.INTEGER(unsigned=True), primary_key=True)

    @declared_attr
    def region_id(cls):
        return Column("regionID", ForeignKey("regions.uid"),
                      nullable=False)

    @declared_attr
    def source_id(cls):
        return Column("sourceID", ForeignKey("sources.uid"),
                      nullable=False)

    # methods
    @classmethod
    def make_query(cls, session, query):
        if (query is not None):
            return query
        q = session.query(cls, Region, Source)\
            .join()\
            .filter(Region.uid == cls.region_id, Source.uid == cls.source_id,)\
            .with_hint(Source, 'USE INDEX (PRIMARY)')\
            .with_hint(Region, 'USE INDEX (PRIMARY)')
        return q

    @classmethod
    def select_all(cls, session, query):
        q = cls.make_query(session, query)
        return q

    @classmethod
    def select_by_chrom(cls, session, query, chrom):
        """
        Query objects by genomic location, 
        retrieve all objects that are in a sepcific chrom.
        """
        q = query.filter(Region.chrom == chrom)

        return q

    @classmethod
    def select_by_overlapping_location(cls, session, query, chrom, start, end):
        """
        Query objects by genomic location, 
        retrieve all objects that overlap with range.
        """
        bins = Region._compute_bins(start, end)
        q = query.filter(Region.chrom == chrom,
                         Region.start < end,
                         Region.end > start)\
            .filter(Region.bin.in_(bins))

        return q

    @classmethod
    def select_by_within_location(cls, session, query, chrom, start, end):
        """
        Query objects by genomic location, 
        retrieve all objects that are within range.
        """
        bins = Region._compute_bins(start, end)
        print(bins)
        q = query.filter(Region.chrom == chrom,
                         Region.start >= start,
                         Region.end <= end)\
            .filter(Region.bin.in_(bins))

        return q

    @classmethod
    def select_by_exact_location(cls, session, query, chrom, start, end):
        """
        Query objects by exact genomic location.
        """
        bins = Region._compute_bins(start, end)
        q = query.filter(Region.bin.in_(bins))\
            .filter(Region.chrom == chrom,
                    Region.start == start,
                    Region.end == end)
        return q

    @classmethod
    def select_by_location(cls, session, query, chrom, start=None, end=None, location=None):
        """
        Query objects by genomic location.
        """
        q = cls.make_query(session, query)
        if (start is None and end is None and location is None): 
            q = cls.select_by_chrom(session, q, chrom)
        elif location == 'exact':
            q = cls.select_by_exact_location(session, q,  chrom, start, end)
        elif location == 'within':
            q = cls.select_by_within_location(session, q, chrom, start, end)
        elif location == 'overlapping':
            q = cls.select_by_overlapping_location(
                session, q,  chrom, start, end)
        else:
            return False
        return q

    @classmethod
    def select_by_uids(cls, session, query, uids):
        """
        filter query by uids.
        """
        q = cls.make_query(session, query)
        q = q.filter(cls.uid.in_(uids))
        return q

    @classmethod
    def select_by_sources(cls, session, query, sources):
        """
        filter query by sources.
        """
        q = cls.make_query(session, query)
        q = q.filter(Source.name.in_(sources))
        return q

    @classmethod
    def is_unique(cls, session, regionID, sourceID):
        """
        Check by unique condition. Minimum condition is is regionID and sourceID.
        Additional conditions require overriding of method in child class.
        """
        q = session.query(cls).filter(cls.region_id == regionID,
                                      cls.source_id == sourceID)
        q = q.all()
        return len(q) == 0

    @classmethod
    def as_genomic_feature(self, feat):
        """
        Return feature as GenomicFeature object. 
        """
        return GenomicFeature(
            feat.Region.chrom,
            int(feat.Region.start),
            int(feat.Region.end),
            feat_type=self.__tablename__,
            feat_id="%s_%s" % (self.__tablename__,
                               getattr(feat, self.__name__).uid),
            qualifiers=None
        )
