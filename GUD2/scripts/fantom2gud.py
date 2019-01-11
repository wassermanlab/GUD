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
from GUD2.ORM.enhancer import Enhancer
from GUD2.ORM.experiment import Experiment
from GUD2.ORM.region import Region
from GUD2.ORM.sample import Sample
from GUD2.ORM.source import Source
from GUD2.ORM.tss import TSS

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
    samples = {}
    sample_ids = []
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
        samples.setdefault(sample_id, {
            "sample_name": sample_name,
            "cell_or_tissue": cell_or_tissue,
            "add": add,
            "treatment": treatment,
            "cell_line": cell_line,
            "cancer": cancer
        })

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

    # Create enhancer/TSS tables
    if feat_type == "enhancer":
        table = Enhancer()
        lines = GUDglobals.parse_csv_file(matrix_file, gz)
        counts_start_at = 1
    if feat_type == "tss":
        table = TSS()
        lines = GUDglobals.parse_tsv_file(matrix_file, gz)
        counts_start_at = 7    
    if not engine.has_table(feat_type):
        table.metadata.bind = engine
        table.metadata.create_all(engine)
    # For each line...
    for line in lines:
        # Skip comments
        if line[0].startswith("#"): continue
        # Get sample IDs
        if len(sample_ids) == 0:
            for sample in line[counts_start_at:]:
                m = re.search("(CNhs\d+)", unquote(sample))
                sample_ids.append(m.group(1))
        # If enhancer/TSS...
        elif line[0].startswith("chr") or line[0].startswith("\"chr"):
            # Initialize
            data = {}
            rows = []
            total_tpm = 0.0
            # Get coordinates
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
                if n:
                    gene = n.group(2)
                    tss_id = n.group(1)
                else:
                    gene = None
                    tss_id = 1
            # Ignore non-standard chroms, scaffolds, etc.
            m = re.search("^chr(\S+)$", chrom)
            if not m.group(1) in GUDglobals.chroms: continue
            # Get region
            region = Region()
            reg = region.select_by_exact_location(session, chrom, start, end)
            if not reg:
                # Insert region
                region.bin = assign_bin(start, end)
                region.chrom = chrom
                region.start = start
                region.end = end
                session.add(region)
                session.commit()
                reg = region.select_by_exact_location(session, chrom, start, end)
            # For each sample...
            for i in range(counts_start_at, len(line)):
                # Skip sample
                sample_id = sample_ids[i - counts_start_at]
                if not samples[sample_id]["add"]: continue
                # Get data
                name = samples[sample_id]["cell_or_tissue"]
                treatment = samples[sample_id]["treatment"]
                cell_line = samples[sample_id]["cell_line"]
                cancer = samples[sample_id]["cancer"]
                data.setdefault((name, treatment, cell_line, cancer), [])
                data[(name, treatment, cell_line, cancer)].append(float(line[i]))
             # For each sample...
            for name, treatment, cell_line, cancer in data:
                # Get total TPMs in "normal" samples
                if not treatment and not cell_line and not cancer:
                    total_tpm += float(sum(data[name, treatment, cell_line, cancer]) /
                        len(data[name, treatment, cell_line, cancer]))
            # For each sample...
            for name, treatment, cell_line, cancer in data:
                # Skip if enhancer/TSS not expressed in sample
                avg_tpm = float(sum(data[name, treatment, cell_line, cancer]) /
                    len(data[name, treatment, cell_line, cancer]))
                if avg_tpm == 0: continue
                # Get sample
                sample = Sample()
                sam = sample.select_by_exact_sample(session, name, treatment, cell_line, cancer)
                if not sam:    
                    sample.name = name
                    sample.treatment = treatment
                    sample.cell_line = cell_line
                    sample.cancer = cancer
                    session.add(sample)
                    session.commit()
                    sam = sample.select_by_exact_sample(session, name, treatment, cell_line, cancer)
                if feat_type == "enhancer":
                    enhancer = Enhancer()
                    if enhancer.is_unique(session, reg.uid, sou.uid, sam.uid, exp.uid):
                        enhancer.regionID = reg.uid
                        enhancer.sourceID = sou.uid
                        enhancer.sampleID = sam.uid
                        enhancer.experimentID = exp.uid
                        rows.append(enhancer)
                if feat_type == "tss":
                    tss = TSS()
                    if tss.is_unique(session, reg.uid, sou.uid, sam.uid, exp.uid):
                        tss.regionID = reg.uid
                        tss.sourceID = sou.uid
                        tss.sampleID = sam.uid
                        tss.experimentID = exp.uid
                        tss.strand = strand
                        tss.gene = gene
                        tss.tss = tss_id
                        tss.avg_tpm = "%.3f" % avg_tpm
                        tss.rel_tpm = "%.3f" % (avg_tpm * 100.0 / total_tpm)
                        rows.append(tss)
            session.add_all(rows)
            session.commit()

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