from sqlalchemy import (Column, ForeignKey)
from sqlalchemy.dialects import mysql
from .genomic_feature import GenomicFeature
from .region import Region
from .source import Source
from .chrom import Chrom
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
    def get_last_uid_region(cls, session, chrom, start, end):
        # returns last uid from region
        if (chrom is None and start is None and end is None):
            return 0

        q = session.query(cls)\
            .join(Region, Region.uid == cls.region_id)\
            .filter(Region.chrom == chrom)

        bins = Region._compute_bins(start, end)
        q = q.filter(Region.bin.in_(bins))
        print(q.order_by(cls.uid).limit(1).statement.compile(compile_kwargs={"literal_binds": True}))
    
        res = q.order_by(cls.uid).limit(1).first()
        if (res == None):
            return None
        else:
            res = res.uid
        print (res-1)
        return (res-1)

    @classmethod
    def get_unique_source_names(cls, session):
        source_ids = session.query(cls.source_id).distinct().all()
        source_ids = [s[0] for s in source_ids]
        source_names = session.query(Source.name).filter(Source.uid.in_(source_ids)).distinct().all()
        source_names = [s[0] for s in source_names]
        return source_names

    @classmethod
    def make_query(cls, session, query):
        if (query is not None):
            return query
        q = session.query(cls, Region.chrom, Region.start, Region.end, Source.name.label('sourceName'))\
            .prefix_with("STRAIGHT_JOIN")\
            .join(Region, Region.uid == cls.region_id)\
            .join(Source, Source.uid == cls.source_id)
        return q

    @classmethod
    def select_all(cls, session, query):
        q = cls.make_query(session, query)
        return q

    @classmethod
    def select_by_overlapping_location(cls, session, query, chrom, start, end):
        """
        Query objects by genomic location, 
        retrieve all objects that overlap with range.
        """
        # bins = Region._compute_bins(start, end)
        # q = query.filter(Region.chrom == chrom,
        #                  Region.start < end,
        #                  Region.end > start)\
        #     .filter(Region.bin.in_(bins))
        regionIDs = cls._get_region_ids_in_location(session, chrom, start, end,
                                                    location="overlapping")
        q = query.filter(cls.region_id.in_(regionIDs))    
       
        return q

    @classmethod
    def select_by_within_location(cls, session, query, chrom, start, end):
        """
        Query objects by genomic location, 
        retrieve all objects that are within range.
        """
        # bins = Region._compute_bins(start, end)
        # q = query.filter(Region.bin.in_(bins))\
        #     .filter(Region.chrom == chrom,
        #                  Region.start >= start,
        #                  Region.end <= end)\
        regionIDs = cls._get_region_ids_in_location(session, chrom, start, end,
                                                    location="within")
        q = query.filter(cls.region_id.in_(regionIDs))    

        return q

    @classmethod
    def select_by_exact_location(cls, session, query, chrom, start, end):
        """
        Query objects by exact genomic location.
        """
        # bins = Region._compute_bins(start, end)
        # q = query.filter(Region.bin.in_(bins))\
        #     .filter(Region.chrom == chrom,
        #             Region.start == start,
        #             Region.end == end)
        regionIDs = cls._get_region_ids_in_location(session, chrom, start, end,
                                                    location="exact")
        q = query.filter(cls.region_id.in_(regionIDs))    
       
        return q   

    @classmethod
    def select_by_location(cls, session, query, chrom, start=None, end=None,
                           location="within"):
        """
        Query objects by genomic location.
        """
        q = cls.make_query(session, query)
        if location == "exact":
            q = cls.select_by_exact_location(session, q,  chrom, start, end)
        elif location == "within":
            q = cls.select_by_within_location(session, q, chrom, start, end)
        elif location == "overlapping":
            q = cls.select_by_overlapping_location(session, q,  chrom, start,
                                                   end)
        else:
            return False

        return q

    @classmethod
    def select_by_chrom(cls, session, query, chrom):
        """
        Query objects by genomic location, 
        retrieve all objects that are within range.
        """
        q = query.filter(Region.chrom == chrom)

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
            feat.chrom,
            int(feat.start),
            int(feat.end),
            feat_type=self.__tablename__,
            feat_id="%s_%s" % (self.__tablename__,
                               getattr(feat, self.__name__).uid),
            qualifiers=None
        )

    @classmethod
    def _get_region_ids_in_location(cls, session, chrom, start, end,
                                    location="within"):
        regions = Region.select_by_bin_range_and_location(session, chrom, start,
                                                          end, bins=[],
                                                          compute_bins=True,
                                                          location=location)
        return [r.uid for r in regions]