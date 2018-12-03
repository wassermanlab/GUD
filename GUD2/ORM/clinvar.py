from sqlalchemy import (
    Column, Index, PrimaryKeyConstraint, String, ForeignKey,
    UniqueConstraint, CheckConstraint, Integer, Float
)
from sqlalchemy.dialects import mysql
from GUD2.ORM.region import Region
from GUD2.ORM.source import Source
from GUD2.ORM.base import Base
from binning import containing_bins, contained_bins, assign_bin

class ClinVar(Base):

    __tablename__ = "clinvar"

    uid = Column("uid", mysql.INTEGER(unsigned=True))
    regionID = Column("regionID", Integer, ForeignKey('regions.uid'), nullable=False)
    sourceID = Column("sourceID", Integer, ForeignKey('sources.uid'), nullable=False)
    ##fields
    ref = Column("ref", String, nullable=False)
    alt = Column("alt", String, nullable=False)
    clinvarID = Column("clinvarID", String, nullable=False)
    ##info
    ANN_Allele = Column("ANN_Allele", String)
    ANN_Annotation = Column("ANN_Annotation", String)
    ANN_Annotation_Impact = Column("ANN_Annotation_Impact", String)
    ANN_Gene_Name = Column("ANN_Gene_Name", String)
    ANN_Gene_ID = Column("ANN_Gene_ID", String)
    ANN_Feature_Type = Column("ANN_Feature_Type", String)
    ANN_Feature_ID = Column("ANN_Feature_ID", String)
    CADD = Column("CADD", Float)
    CLNDISDB = Column("CLNDISDB", String)
    CLNDN = Column("CLNDN", String)
    CLNSIG = Column("CLNSIG", String)
    gnomad_exome_af_global = Column("gnomad_exome_af_global", Float)
    gnomad_exome_hom_global = Column("gnomad_exome_hom_global", Float)
    gnomad_genome_af_global = Column("gnomad_genome_af_global", Float)
    gnomad_genome_hom_global = Column("gnomad_genome_hom_global", Float) 

    __table_args__ = (
        PrimaryKeyConstraint(uid),
        UniqueConstraint(clinvarID),

        Index("ix_clinvar", regionID),
        Index("ix_clinvar_clnsig", CLNSIG),
        Index("ix_clinvar_feature", ANN_Feature_Type),
        Index("ix_clinvar_gene_name", ANN_Gene_Name),
        Index("ix_clinvar_id", clinvarID),

        {
            "mysql_engine": "MyISAM",
            "mysql_charset": "utf8"
        }
    )

    @classmethod
    def select_by_location(cls, session, chrom, start, end):
        """
        Query objects based off of their location being within the start only
        motifs through that  
         """
        # print(chrom, start, end)
        # q = Region.select_by_bin_range(session, chrom, start, end, [], True, False)
        bin = assign_bin(start, end)
        q = session.query(cls, Region).\
        join().\
        filter(Region.uid == cls.regionID).\
        filter(Region.chrom == chrom, Region.end > start, Region.start < end).\
        filter(Region.bin == bin)
        return q.all() 

    @classmethod
    def is_unique(cls, session, clinvarID):
        q = session.query(cls).filter(cls.clinvarID == clinvarID)
        return len(q.all()) == 0

    # def __str__(self):
    #     return "{}\t{}".format(self.motif, self.pathogenicity)

    # def __repr__(self):
    #     return "<ShortTandemRepeat(uid={}, regionID={}, sourceID={}, motif={}, pathogencity={})>".format(
    #         self.uid, self.regionID, self.sourceID, self.motif, self.pathogenicity)
