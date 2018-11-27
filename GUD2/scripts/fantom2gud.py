#!/usr/bin/env python

import os, sys, re
import argparse
from binning import assign_bin
import getpass
from sqlalchemy import create_engine
from sqlalchemy.orm import mapper, scoped_session, sessionmaker
from sqlalchemy_utils import database_exists
from urllib2 import unquote
import warnings

# Import from GUD module
from GUD2 import GUDglobals
#from GUD2.ORM.enhancer import Enhancer ##todo
from GUD2.ORM.tss import TSS
from GUD2.ORM.region import Region 
from GUD2.ORM.source import Source 
from GUD2.ORM.sample import Sample 
from GUD2.ORM.experiment import Experiment 
from GUD2.scripts.sample_names import sample_names
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

    feats = ["enhancer", "tss"]
    parser.add_argument("feat_type", choices=feats, help="Type of genomic feature", metavar="feature_type")

    # Optional args
    parser.add_argument("-k", "--keep", action="store_true", help="Keep original sample fanmes from FANTOM (default = False)")
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

    # Set default
    if not args.db:
        args.db = args.genome

    return args

def insert_fantom_to_gud_db(user, host, port, db, matrix_file,
    feat_type, source_name, keep=False):

    # Initialize
    coordinates = set()
    fantom_sample_names = []
    db_name = "mysql://{}:@{}:{}/{}".format(
        user, host, port, db)
    if not database_exists(db_name):
        raise ValueError("GUD db does not exist: %s" % db_name)
    session = scoped_session(sessionmaker())
    engine = create_engine(db_name, echo=False)
    session.remove()
    session.configure(bind=engine, autoflush=False,
        expire_on_commit=False)
    if matrix_file.endswith(".gz"): gz = True
    else: gz = False

    # Initialize table
    if feat_type == "enhancer":
        table = Enhancer()
    if feat_type == "tss":
        table = TSS()
    if not engine.has_table(table.__tablename__):
        try:
            table.metadata.bind = engine
            table.metadata.create_all(engine)
        except:
            raise ValueError("Cannot create table: %s" % table.__tablename__)

    # For each line...
    if feat_type == "enhancer":
        lines = GUDglobals.parse_csv_file(matrix_file, gz)
        counts_start_at = 1
    if feat_type == "tss":
        lines = GUDglobals.parse_tsv_file(matrix_file, gz)
        counts_start_at = 7
    for line in lines:
        # Skip comments
        if line[0].startswith("#"): continue
        # If no samples...
        if len(fantom_sample_names) == 0:
            for sample in line[counts_start_at:]:
                fantom_sample_names.append(unquote(sample))
        # ... Else...
        elif line[0].startswith("chr") or line[0].startswith("\"chr"):
            # Initialize
            samples = {}
            total_tpms = 0.0
            # Get chrom, start, end
            if feat_type == "enhancer":
                feat = Enhancer()
                m = re.search("(chr\S+)\:(\d+)\-(\d+)", line[0])
            if feat_type == "tss":
                feat = TSS()
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
            ## Region
            region = Region()
            reg = region.select_by_pos(session, chrom, start, end)
            if not reg: 
                region.bin = assign_bin(int(start), int(end))
                region.chrom = chrom
                region.start = start
                region.end = end 
                session.merge(region)
                session.commit()
                reg = region.select_by_pos(session, chrom, start, end)
            ## Source
            source = Source()
            sou = source.select_by_name(session, source_name)
            if not sou: 
                source.name = source_name
                session.merge(source)
                session.commit()
                sou = source.select_by_name(session, source_name)
            ## Experiment
            experiment = Experiment()
            ex = experiment.select_by_name(session, "CAGE")
            if not ex:
                experiment.name = "CAGE"
                session.merge(experiment)
                session.commit()
                ex = experiment.select_by_name(session, "CAGE")

            if feat_type == "tss":
                feat.gene = "%s:%s-%s,%s" % (chrom, start + 1, end, strand)
                feat.tss = 1
                feat.strand = strand
                if n:
                    feat.gene = n.group(2)
                    feat.tss = n.group(1)
            # For each sample...
            for i in range(counts_start_at, len(line)):
                # Initialize
                tpm = "%.3f" % float(line[i])
                sample = fantom_sample_names[i - counts_start_at]
                # Keep original sample names
                if keep:
                    samples.setdefault(sample, [float(tpm)])
                    total_tpms += float(tpm)
                else:
                    m = re.search("(CNhs\d+)", sample)
                    if m.group(1) in sample_names:
                        samples.setdefault(sample_names[m.group(1)], [])
                        samples[sample_names[m.group(1)]].append(float(tpm))
                        total_tpms += float(tpm)
            # For each sample...
            for sample in samples:
                ## Region
                sample = Sample()
                samp = sample.select_exact_sample(session, name, treatment, cell_line, cancer)
                if not reg: 
                    samp.name = 
                    samp.treatment = 
                    samp.cell_line = 
                    samp.cancer = 
                    session.merge(samp)
                    session.commit()
                    samp = sample.select_exact_sample(session, name, treatment, cell_line, cancer)
                model.cell_or_tissue = sample ## ask oriol about this 
                # Skip enhancers with 0 tpm
                if sum(samples[sample]) == 0: continue
                if feat_type == "tss":
                    feat.avg_tpm = "%.3f" % float(sum(samples[sample]) / len(samples[sample]))
                
                ## Sample 
                ## Feature commiting 
                feat.regionID = reg[0].uid
                feat.sourceID = sou[0].uid
                feat.sampleID = samp[0].uid 
                feat.experimentID = ex[0].uid
                session.merge(feat)
                session.commit()

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    # Parse arguments
    args = parse_args()

    # Insert FANTOM data to GUD database
    insert_fantom_to_gud_db(args.user, args.host, args.port,
        args.db, args.matrix, args.feat_type, args.source,
        args.keep)