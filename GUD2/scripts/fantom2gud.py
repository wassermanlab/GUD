#!/usr/bin/env python

import os, sys, re
import argparse
from binning import assign_bin
from datetime import date
import getpass
import pybedtools
from sqlalchemy import create_engine
from sqlalchemy.orm import mapper, scoped_session, sessionmaker
from sqlalchemy_utils import create_database, database_exists
try: from urllib2 import unquote
except: from urllib.parse import unquote
import warnings

# Import from GUD module
from GUD2 import GUDglobals
from GUD2.ORM.tss import TSS
from GUD2.ORM.experiment import Experiment
from GUD2.ORM.region import Region
from GUD2.ORM.sample import Sample
from GUD2.ORM.source import Source

#-------------#
# Functions   #
#-------------#

def parse_args():
    """
    This function parses arguments provided via the command
    line using argparse.
    """

    parser = argparse.ArgumentParser(description="this script inserts \"enhancer\" or \"tss\" data from FANTOM into GUD. download \"hg19_permissive_enhancers_expression_rle_tpm.csv.gz\" and \"hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt.gz\" for enhancer and tss data, respectively.")

    parser.add_argument("matrix", help="Expression (TPM/RLE normalized) matrix across all FANTOM libraries")
    parser.add_argument("samples", help="FANTOM samples (manually-curated)")

    feats = ["enhancer", "tss"]
    parser.add_argument("feat_type", choices=feats, help="Type of genomic feature", metavar="feature_type")

    # Optional args
    parser.add_argument("--source", default="FANTOM", help="Source name (e.g. \"PMID:24670763\" for TSSs or \"PMID:24670764\" for enhancers; default = \"FANTOM\")")

    # MySQL args
    mysql_group = parser.add_argument_group("mysql arguments")
    mysql_group.add_argument("-d", "--db", default="hg19",
        help="Database name (default = \"hg19\")")
    mysql_group.add_argument("-H", "--host", default="localhost",
        help="Host name (default = localhost)")
    mysql_group.add_argument("-P", "--port", default=5506, type=int,
        help="Port number (default = 5506)")
    mysql_group.add_argument("-u", "--user", default=getpass.getuser(),
        help="User name (default = current user)")

    args = parser.parse_args()

    return args

def insert_fantom_to_gud_db(user, host, port, db,
    matrix_file, samples_file, feat_type, source_name):

    # Initialize
    sample_ids = []
    sample_metadata = {}
    db_name = "mysql://{}:@{}:{}/{}".format(
        user, host, port, db)
    if not database_exists(db_name):
        create_database(db_name)
    session = scoped_session(sessionmaker())
    engine = create_engine(db_name, echo=False)
    session.remove()
    session.configure(bind=engine, autoflush=False,
        expire_on_commit=False)
    if matrix_file.endswith(".gz"): gz = True
    else: gz = False

    # Get samples
    for line in GUDglobals.parse_tsv_file(samples_file):
        # Initialize 
        sample_id = line[0]
        sample_name = line[1]
        cell_or_tissue = line[2]
        if line[3] == "Yes": add = True
        else: add = False
        if line[4] == "Yes": treatment = True
        else: treatment = False
        if line[5] == "Yes": cell_line = True
        else: cell_line = False
        if line[6] == "Yes": cancer = True
        else: cancer = False
        # Get sample
        sample_metadata.setdefault(sample_id, {
            "sample_name": sample_name,
            "cell_or_tissue": cell_or_tissue,
            "add": add,
            "treatment": treatment,
            "cell_line": cell_line,
            "cancer": cancer
        })

    # Initialize enhancer table
    if feat_type == "enhancer":
        table = Enhancer()
        lines = GUDglobals.parse_csv_file(matrix_file, gz)
        counts_start_at = 1
    # Initialize tss table
    if feat_type == "tss":
        table = TSS()
        lines = GUDglobals.parse_tsv_file(matrix_file, gz)
        counts_start_at = 7
    # Create table    
    if not engine.has_table(feat_type):
        table.metadata.bind = engine
        table.metadata.create_all(engine)

    # Get experiment
    experiment = Experiment()
    exp = experiment.select_by_name(session, "CAGE")
    if not exp:    
        experiment.name = "CAGE"
        session.add(experiment)
    session.commit()
    exp = experiment.select_by_name(session, "CAGE")

    # Get source
    source = Source()
    sou = source.select_by_name(session, source_name)
    if not sou:    
        source.name = source_name
        session.add(source)
    session.commit()
    sou = source.select_by_name(session, source_name)

    print(exp, sou)
    exit(0)


    # Insert data
    for line in lines:
        # Skip comments
        if line[0].startswith("#"): continue
        # If no samples...
        if len(sample_ids) == 0:
            for sample in line[counts_start_at:]:
                m = re.search("(CNhs\d+)", unquote(sample))
                sample_ids.append(m.group(1))
        # If enhancer/TSS...
        elif line[0].startswith("chr") or line[0].startswith("\"chr"):
            # Initialize
            samples = {}
            total_tpms = 0.0
            # Get chrom, start, end
            if feat_type == "enhancer":
                m = re.search("(chr\S+)\:(\d+)\-(\d+)", line[0])
            if feat_type == "tss":
                m = re.search("(chr\S+)\:(\d+)\.\.(\d+),(\S)", line[0])
                n = re.search("p(\d+)@(\w+)", line[1])
            chrom = m.group(1)
            start = int(m.group(2))
            end = int(m.group(3))
            if feat_type == "tss":
                strand = m.group(4)
            # Ignore non-standard chroms, scaffolds, etc.
            m = re.search("^chr(\S+)$", chrom)
            if not m.group(1) in GUDglobals.chroms: continue

            # Get region
            region = Region()
            reg = region.select_by_pos(session, chrom, start, end)
            if not reg:
                region.bin = assign_bin(start, end)
                region.chrom = chrom
                region.start = start
                region.end = end
                session.add(region)
            reg = region.select_by_pos(session, chrom, start, end)

            if feat_type == "enhancer":
                feature = Enhancer()
            if feat_type == "tss":
                feature = TSS()
            if tss.is_unique(session, reg.uid, sou.uid):
                tss.score = line[6]
                tss.regionID = reg.uid
                tss.sourceID = sou.uid
                rows.append(conservation)

                tss.regionID = Column("regionID", Integer, ForeignKey('regions.uid'), nullable=False)
                tss.sourceID = Column("sourceID", Integer, ForeignKey('sources.uid'), nullable=False)
                tss.sampleID = Column("sampleID", Integer, ForeignKey('samples.uid'), nullable=False)
                tss.experimentID = Column("experimentID", Integer, ForeignKey('experiments.uid'), nullable=False)
                tss.gene = Column("gene", String(75), ForeignKey('genes.name2'))
                tss.tss = Column("tss", mysql.INTEGER(unsigned=True))
                tss.avg_tpm = Column("avg_tpm", Float, nullable=False)
                tss.strand = Column("strand", mysql.CHAR(1), nullable=False)
#
#                
#            if len(rows) == 100000:
#                session.add_all(rows)
#                session.commit()
#                rows = []  
#        session.add_all(rows)
#        session.commit()
#
#            #region entry 
#            chrom = line[1]
#            start = int(line[2])
#            end = int(line[3])
#            region = Region()
#            source = Source()
#            reg = region.select_by_pos(session, chrom, start, end)
#            sou = source.select_by_name(session, source_name.group(1))
#
#
#
#            # Create model
#            model = Model()
#            model.bin = assign_bin(int(start), int(end))
#            model.chrom = chrom
#            model.start = start
#            model.end = end
#            model.experiment_type = "CAGE"
#            model.source_name = source_name
#            model.date = today
#            if feat_type == "tss":
#                model.gene = "%s:%s-%s,%s" % (chrom, start + 1, end, strand)
#                model.tss = 1
#                model.strand = strand
#                if n:
#                    model.gene = n.group(2)
#                    model.tss = n.group(1)
#            # For each sample...
#            for i in range(counts_start_at, len(line)):
#                # Initialize
#                tpm = "%.3f" % float(line[i])
#                sample = fantom_sample_names[i - counts_start_at]
#                # Keep original sample names
#                if keep:
#                    samples.setdefault(sample, [float(tpm)])
#                    total_tpms += float(tpm)
#                else:
#                    m = re.search("(CNhs\d+)", sample)
#                    if m.group(1) in sample_names:
#                        samples.setdefault(sample_names[m.group(1)], [])
#                        samples[sample_names[m.group(1)]].append(float(tpm))
#                        total_tpms += float(tpm)
#            # For each sample...
#            for sample in samples:
#                model.cell_or_tissue = sample
#                # Skip enhancers with 0 tpm
#                if sum(samples[sample]) == 0: continue
#                if feat_type == "tss":
#                    model.avg_tpm = "%.3f" % float(sum(samples[sample]) / len(samples[sample]))
##                    model.percent_tpm = "%.3f" % float(sum(samples[sample]) * 100 / total_tpms)
#                # Upsert model & commit
#                session.merge(model)
#                session.commit()

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    # Parse arguments
    args = parse_args()

    # Insert FANTOM data to GUD database
    insert_fantom_to_gud_db(args.user, args.host, args.port,
        args.db, args.matrix, args.samples, args.feat_type,
        args.source)