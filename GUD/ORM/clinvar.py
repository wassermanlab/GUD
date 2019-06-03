from binning import (
    containing_bins,
    contained_bins
)
from sqlalchemy import (
    Column, Index, PrimaryKeyConstraint, String, ForeignKey,
    UniqueConstraint, CheckConstraint, Integer, Float
)
from sqlalchemy.dialects import mysql

from .base import Base
from .region import Region
from .source import Source

class ClinVar(Base):

    __tablename__ = "clinvar"

    uid = Column("uid", mysql.INTEGER(unsigned=True))
    regionID = Column("regionID", Integer, ForeignKey('regions.uid'), nullable=False)
    sourceID = Column("sourceID", Integer, ForeignKey('sources.uid'), nullable=False)
    ##fields
    ref = Column("ref", mysql.VARCHAR(5000), nullable=False) ## check these
    alt = Column("alt", mysql.VARCHAR(5000), nullable=False) ## check these
    clinvarID = Column("clinvarID", mysql.VARCHAR(7), nullable=False)
    ##info
    ANN_Annotation = Column("ANN_Annotation", mysql.VARCHAR(250))
    ANN_Annotation_Impact = Column("ANN_Annotation_Impact", mysql.VARCHAR(50))
    ANN_Gene_Name = Column("ANN_Gene_Name", mysql.VARCHAR(50))
    ANN_Gene_ID = Column("ANN_Gene_ID", mysql.VARCHAR(50))
    ANN_Feature_Type = Column("ANN_Feature_Type", mysql.VARCHAR(500))
    ANN_Feature_ID = Column("ANN_Feature_ID", mysql.VARCHAR(50))
    CADD = Column("CADD", Float)
    CLNDISDB = Column("CLNDISDB", mysql.VARCHAR(3000)) ## check 
    CLNDN = Column("CLNDN", mysql.VARCHAR(3000))
    CLNSIG = Column("CLNSIG", mysql.VARCHAR(500))
    gnomad_exome_af_global = Column("gnomad_exome_af_global", Float)
    gnomad_exome_hom_global = Column("gnomad_exome_hom_global", Float)
    gnomad_genome_af_global = Column("gnomad_genome_af_global", Float)
    gnomad_genome_hom_global = Column("gnomad_genome_hom_global", Float) 

    __table_args__ = (
        PrimaryKeyConstraint(uid),
        UniqueConstraint(clinvarID),

        Index("ix_clinvar", regionID),
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
        bins = set(containing_bins(start, end) + contained_bins(start, end))
        q = session.query(cls, Region).\
        join().\
        filter(Region.uid == cls.regionID).\
        filter(Region.chrom == chrom, Region.end > start, Region.start < end).\
        filter(Region.bin.in_(bins))
        return q.all()

    @classmethod
    def select_by_name(cls, session, clinvarID):
        """
        Query refGene objects by common name. If no name is provided,
        query all genes.
        """
        q = session.query(cls)
        q = session.query(cls, Region).\
        join().\
        filter(Region.uid == cls.regionID).\
        filter(cls.clinvarID == clinvarID)
        return q.first()

    @classmethod
    def is_unique(cls, session, clinvarID):
        q = session.query(cls).filter(cls.clinvarID == clinvarID)
        return len(q.all()) == 0

    def __str__(self):
        return "REF\tALT\tclinvarID\tannotation\tannotation_impact\tfeature_type\tCLNDISDB\tCLINSIG\n{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".\
        format(self.ref, self.alt, self.clinvarID, 
        self.ANN_Annotation, self.ANN_Annotation_Impact, self.ANN_Feature_Type,
        self.CLNDISDB, self.CLNSIG)

    def __repr__(self):
        return "<ShortTandemRepeat(uid={}, regionID={}, sourceID={},\
        ref={}, alt={}, clinvarID={},\
        ANN_Annotation={}, ANN_Annotation_Impact={}, ANN_Gene_Name={},\
        ANN_Gene_ID={}, ANN_Feature_Type={}, ANN_Feature_ID={},\
        CADD={}, CLNDISDB={}, CLNDN={}, CLNSIG={},\
        gnomad_exome_af_global={}, gnomad_exome_hom_global={}, gnomad_genome_af_global={}, gnomad_genome_hom_global={})>".format(
            self.uid, self.regionID, self.sourceID, self.ref, self.alt, self.clinvarID,
            self.ANN_Annotation, self.ANN_Annotation_Impact, self.ANN_Gene_Name,
            self.ANN_Gene_ID, self.ANN_Feature_Type, self.ANN_Feature_ID,
            self.CADD, self.CLNDISDB, self.CLNDN, self.CLNSIG, 
            self.gnomad_exome_af_global, self.gnomad_exome_hom_global,
            self.gnomad_genome_af_global, self.gnomad_genome_hom_global)
