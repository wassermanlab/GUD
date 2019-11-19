from sqlalchemy import (Column, Index, String, ForeignKey, UniqueConstraint) 
import re 
from sqlalchemy.dialects import mysql
from .base import Base
from .experiment import Experiment
from .gene import Gene
from .region import Region
from .source import Source
from .genomicFeatureMixin2 import GFMixin2
from sqlalchemy.ext.declarative import declared_attr


class TSS(GFMixin2, Base):
    # table declerations
    __tablename__ = "transcription_start_sites"

    gene = Column("gene", String(75), ForeignKey("genes.name2"))
    tss = Column("tss", mysql.INTEGER(unsigned=True))

    @declared_attr
    def sample_id(cls):
        return Column("sampleIDs", mysql.LONGBLOB, nullable=False)

    avg_expression_levels = Column("avg_expression_levels", mysql.LONGBLOB,
                                   nullable=False)

    @declared_attr
    def __table_args__(cls):
        return (
            UniqueConstraint(cls.region_id, cls.sample_id, cls.experiment_id, cls.gene, cls.tss),
            Index("ix_join", cls.region_id, cls.experiment_id, cls.source_id),
            Index("ix_gene_tss", cls.gene, cls.tss),
            {"mysql_engine": "InnoDB", "mysql_charset": "utf8"}
        )

    # class methods 
    @classmethod
    def is_unique(cls, session, regionID, sourceID,
                  experimentID, gene, tss):

        q = session.query(cls)\
            .filter(cls.region_id == regionID, cls.source_id == sourceID,
                    cls.experiment_id == experimentID, 
                    cls.gene == gene, cls.tss == tss)

        return len(q.all()) == 0

    @classmethod
    def select_by_genes(cls, session, query, genes):
        """
        select by genes
        """
        q = cls.make_query(session, query)
        q = q.filter(cls.gene.in_(genes))

        return q

    @classmethod 
    def select_all_genic_tss(cls, session, query):
        """
        Query all objects associated with a gene.
        """
        q = cls.make_query(session, query)
        q = q.filter(cls.gene != None)

        return q

    @classmethod
    def as_genomic_feature(self, feat):

        # Initialize
        isfloat = re.compile("\d+(\.\d+)?")
        sampleIDs = []
        avg_expression_levels = []

        # For each exon start...
        for i in str(feat.TSS.sample_id).split(","):
            if i.isdigit():
                sampleIDs.append(int(i))

        # For each exon end...
        for i in str(feat.TSS.avg_expression_levels).split(","):
            if isfloat.match(i):
                avg_expression_levels.append(float(i))

        # Define qualifiers
        try:
            experiment_name = feat.Experiment.name
        except:
            experiment_name = "CAGE"
        qualifiers = {
            "uid": feat.TSS.uid,
            "gene": feat.TSS.gene,
            "tss": feat.TSS.tss,
            "sampleIDs": feat.TSS.sample_id,
            "avg_expression_levels": feat.TSS.avg_expression_levels,
            "experiment": experiment_name,
            "source": feat.Source.name,
        }

        genomic_feature = super().as_genomic_feature(feat)
        genomic_feature.qualifiers = qualifiers
        return genomic_feature
