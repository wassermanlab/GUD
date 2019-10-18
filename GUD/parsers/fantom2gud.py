#!/usr/bin/env python

import argparse
from binning import assign_bin
from copy import copy
from functools import partial
import getpass
from multiprocessing import cpu_count
import os
import re
import subprocess
import sys
import tarfile
# Python 3+
if sys.version_info > (3, 0):
    from urllib.request import urlretrieve
    from urllib.parse import unquote
# Python 2.7
else:
    from urllib import urlretrieve
    from urllib2 import unquote

# Import from GUD module
from GUD import GUDUtils
from GUD.ORM.enhancer import Enhancer
from GUD.ORM.experiment import Experiment
from GUD.ORM.expression import Expression
from GUD.ORM.region import Region
from GUD.ORM.sample import Sample
from GUD.ORM.source import Source
from GUD.ORM.tss import TSS
from . import ParseUtils

usage_msg = """
usage: %s --genome STR --samples FILE --feature STR [-h] [options]
""" % os.path.basename(__file__)

help_msg = """%s
inserts features from FANTOM5 (Functional Annotation of the
Mammalian Genome) into GUD.

  --genome STR        genome assembly
  --samples FILE      FANTOM5 samples (manually-curated)
  --feature STR       type of genomic feature ("enhancer"
                      "tss")

optional arguments:
  -h, --help          show this help message and exit
  --dummy-dir DIR     dummy directory (default = "/tmp/")
  -t, --test          limit the total of inserts to ~1K per
                      thread for testing (default = False)
  --threads INT       number of threads to use (default = %s)

mysql arguments:
  -d STR, --db STR    database name (default = "%s")
  -H STR, --host STR  host name (default = "localhost")
  -p STR, --pwd STR   password (default = ignore this option)
  -P STR, --port STR  port number (default = %s)
  -u STR, --user STR  user name (default = current user)
""" % (usage_msg, (cpu_count() - 1), GUDUtils.db, GUDUtils.port)

#-------------#
# Functions   #
#-------------#

def parse_args():
    """
    This function parses arguments provided via the command line and returns an {argparse} object.
    """

    parser = argparse.ArgumentParser(add_help=False)

    # Mandatory args
    parser.add_argument("--genome")
    parser.add_argument("--samples")

    # Optional args
    optional_group = parser.add_argument_group("optional arguments")
    optional_group.add_argument("-h", "--help", action="store_true")
    optional_group.add_argument("--dummy-dir", default="/tmp/")
    optional_group.add_argument("-m", "--merge", action="store_true")
    optional_group.add_argument("-t", "--test", action="store_true")
    optional_group.add_argument("--threads", default=(cpu_count() - 1))

    # MySQL args
    mysql_group = parser.add_argument_group("mysql arguments")
    mysql_group.add_argument("-d", "--db", default=GUDUtils.db)
    mysql_group.add_argument("-H", "--host", default="localhost")
    mysql_group.add_argument("-p", "--pwd")
    mysql_group.add_argument("-P", "--port", default=GUDUtils.port)
    mysql_group.add_argument("-u", "--user", default=getpass.getuser())

    args = parser.parse_args()

    check_args(args)

    return(args)

def check_args(args):
    """
    This function checks an {argparse} object.
    """

    # Print help
    if args.help:
        print(help_msg)
        exit(0)

    # Check mandatory arguments
    if not args.genome or not args.samples or not args.feature:
        error = ["%s\n%s" % (usage_msg, os.path.basename(__file__)), "error", "argument \"--genome\" \"--samples\" \"--feature\" are required\n"]
        print(": ".join(error))
        exit(0)

    # Check "-t" argument
    try:
        args.threads = int(args.threads)
    except:
        error = ["%s\n%s" % (usage_msg, os.path.basename(__file__)), "error", "argument \"-t\" \"--threads\"", "invalid int value", "\"%s\"\n" % args.threads]
        print(": ".join(error))
        exit(0)

    # Check MySQL password
    if not args.pwd:
        args.pwd = ""

def main():

    # Parse arguments
    args = parse_args()

    # Set MySQL options
    GUDUtils.user = args.user
    GUDUtils.pwd  = args.pwd
    GUDUtils.host = args.host
    GUDUtils.port = args.port
    GUDUtils.db = args.db

    # Insert FANTOM5 data
    fantom_to_gud(args.genome, args.samples, args.feature, args.dummy_dir, args.test, args.threads)

def fantom_to_gud(genome, samples_file, feat_type, dummy_dir="/tmp/", test=False, threads=1):

    # Initialize
    global chroms
    global engine
    global experiment
    global Feature
    global Session
    global samples
    global source
    global tpms_start_at
    experiment_type = "CAGE"
    source_name = "FANTOM5"

    # Testing
    if test:
        global current_process
        from multiprocessing import current_process

    # Get database name
    db_name = GUDUtils._get_db_name()

    # Get engine/session
    engine, Session = GUDUtils.get_engine_session(db_name)

    # Initialize parser utilities
    ParseUtils.genome = genome
    ParseUtils.dbname = db_name
    ParseUtils.engine = engine

    # If database does not exist...
    ParseUtils.initialize_gud_db()

    # Create table
    if feat_type == "enhancer":
        Feature = Enhancer
        tpms_start_at = 1
    else:
        Feature = TSS
        tpms_start_at = 7
        ParseUtils.create_table(Expression)
    ParseUtils.create_table(Feature)

    # Start a new session
    session = Session()

    # Get valid chromosomes
    chroms = ParseUtils.get_chroms(session)

    # Get experiment
    experiment = Experiment()
    experiment.name = experiment_type
    ParseUtils.upsert_experiment(session, experiment)
    experiment = ParseUtils.get_experiment(session, experiment_type)    

    # Download data
    download_file = _download_data(genome, dummy_dir)

    # Get samples
    samples = _get_samples(session, samples_file)
    print(samples)

    # Get source
    source = Source()
    source.name = source_name
    ParseUtils.upsert_source(session, source)
    sources = ParseUtils.get_source(session, source_name)
    source = next(iter(sources))

    # This is ABSOLUTELY necessary to prevent MySQL from crashing!
    session.close()
    engine.dispose()

def _download_data(genome, feat_type, dummy_dir="/tmp/"):

    # Python 3+
    if sys.version_info > (3, 0):
        from urllib.request import urlretrieve
    # Python 2.7
    else:
        from urllib import urlretrieve

    # Initialize
    url = "http://fantom.gsc.riken.jp/5/datafiles/"

    # If latest human/mouse genomes...
    if genome == "hg38" or genome == "mm10"
        url += "reprocessed/%s_latest/extra" % genome
        # Get experiment
        if feat_type == "enhancer":
            ftp_file = "F5.%s.enhancers.expression.usage.matrix.gz" %s genome
        else:
            
    else:
        url += "latest/extra"



human_permissive_enhancers_phase_1_and_2_expression_tpm_matrix.txt.gz

    dummy_file = os.path.join(dummy_dir, ftp_file)
    if not os.path.exists(dummy_file):
        f = urlretrieve(os.path.join(url, ftp_file), dummy_file)

    return(dummy_file)

def _get_samples(session, file_name):

    # Initialize
    samples = {}

    # For each line...
    for line in ParseUtils.parse_tsv_file(file_name):

        # If add...
        if line[3] == "Yes":

            # Get sample
            sample = Sample()
            sample.name = line[2]
            sample.treatment = False
            if line[4] == "Yes":
                sample.treatment = True
            sample.cell_line = False
            if line[5] == "Yes":
                sample.cell_line = True
            sample.cancer = False
            if line[6] == "Yes":
                sample.cancer = True

            # Upsert sample
            ParseUtils.upsert_sample(session, sample)

            # Get sample ID
            sample = ParseUtils.get_sample(session, sample.name, sample.X, sample.Y, sample.treatment, sample.cell_line, sample.cancer)

            # Add sample
            samples.setdefault(line[0], sample.uid)

    return(samples)


    # # Get samples
    # for line in GUDglobals.parse_tsv_file(
    #     samples_file
    # ):
    #     # Initialize 
    #     sample_id = line[0]
    #     sample_name = line[1]
    #     cell_or_tissue = line[2]
    #     if line[3] == "Yes": add = True
    #     else: add = False
    #     if line[4] == "Yes": treatment = True
    #     else: treatment = False
    #     if line[5] == "Yes": cell_line = True
    #     else: cell_line = False
    #     if line[6] == "Yes": cancer = True
    #     else: cancer = False
    #     # Get sample
    #     samples.setdefault(sample_id,
    #         {
    #             "sample_name": sample_name,
    #             "cell_or_tissue": cell_or_tissue,
    #             "add": add,
    #             "treatment": treatment,
    #             "cell_line": cell_line,
    #             "cancer": cancer
    #         }
    #     )

    # # Get experiment
    # experiment = Experiment()
    # if experiment.is_unique(session, "CAGE"):    
    #     experiment.name = "CAGE"
    #     session.add(experiment)
    #     session.commit()
    # exp = experiment.select_by_name(
    #     session,
    #     "CAGE"
    # )

    # # Get source
    # source = Source()
    # if source.is_unique(session, source_name):    
    #     source.name = source_name
    #     session.add(source)
    #     session.commit()
    # sou = source.select_by_name(
    #     session,
    #     source_name
    # )

    # # Create table
    # if feat_type == "enhancer":
    #     tpms_start_at = 1
    #     table = Enhancer()
    #     if not engine.has_table(
    #         table.__tablename__
    #     ):
    #         # Create table
    #         table.__table__.create(
    #             bind=engine
    #         )
    #     lines = GUDglobals.parse_csv_file(
    #         matrix_file
    #     )

    # if feat_type == "tss":
    #     tpms_start_at = 7
    #     table = TSS()
    #     if not engine.has_table(
    #         table.__tablename__
    #     ):
    #         # Create table
    #         table.__table__.create(
    #             bind=engine
    #         )
    #     table = Expression()
    #     if not engine.has_table(
    #         table.__tablename__
    #     ):
    #         # Create table
    #         table.__table__.create(
    #             bind=engine
    #         )
    #     lines = GUDglobals.parse_tsv_file(
    #         matrix_file
    #     )

    # # Get valid chromosomes
    # chroms = Chrom.chrom_sizes(session)

    # # For each line...
    # for line in lines:
    #     # Skip comments
    #     if line[0].startswith("#"): continue
    #     # Get sample IDs
    #     if len(sample_ids) == 0:
    #         for sample in line[tpms_start_at:]:
    #             m = re.search(
    #                 "(CNhs\d+)",
    #                 unquote(sample)
    #             )
    #             sample_ids.append(m.group(1))
    #     # If enhancer/TSS...
    #     elif line[0].startswith("chr") or\
    #          line[0].startswith("\"chr"):
    #         # Initialize
    #         data = {}
    #         features = []
    #         sampleIDs = []
    #         avg_expression_levels = []
    #         # Get coordinates
    #         if feat_type == "enhancer":
    #             m = re.search(
    #                 "(chr\S+)\:(\d+)\-(\d+)",
    #                 line[0]
    #             )
    #             chrom = m.group(1)
    #             start = int(m.group(2))
    #             end = int(m.group(3))
    #             strand = None
    #         if feat_type == "tss":
    #             # Initialize
    #             peak_ids = set()
    #             m = re.search(
    #                 "(chr\S+)\:(\d+)\.\.(\d+),(\S)",
    #                 line[0]
    #             )
    #             chrom = m.group(1)
    #             start = int(m.group(2))
    #             end = int(m.group(3))
    #             strand = m.group(4)
    #             for peak in line[1].split(","):
    #                 m = re.search(
    #                     "p(\d+)@(\w+)",
    #                     peak
    #                 )
    #                 if m:
    #                     peak_ids.add((
    #                         m.group(2), m.group(1)
    #                     ))
    #                 else:
    #                     peak_ids.add((None, 1))
    #         # Ignore non-standard chroms,
    #         # scaffolds, etc.
    #         if chrom not in chroms:
    #             continue
    #         # Get region
    #         region = Region()
    #         if region.is_unique(
    #             session,
    #             chrom,
    #             start,
    #             end,
    #             strand
    #         ):
    #             # Insert region
    #             region.bin = assign_bin(start, end)
    #             region.chrom = chrom
    #             region.start = start
    #             region.end = end
    #             region.strand = strand
    #             session.add(region)
    #             session.commit()
    #         reg = region.select_unique(
    #             session,
    #             chrom,
    #             start,
    #             end,
    #             strand
    #         )
    #         # For each sample...
    #         for i in range(tpms_start_at, len(line)):
    #             # Skip sample
    #             sample_id = sample_ids[
    #                 i - tpms_start_at
    #             ]
    #             if not samples[sample_id]["add"]:
    #                 continue
    #             # Get data
    #             name =\
    #                 samples[sample_id]["cell_or_tissue"]
    #             treatment =\
    #                 samples[sample_id]["treatment"]
    #             cell_line =\
    #                 samples[sample_id]["cell_line"]
    #             cancer =\
    #                 samples[sample_id]["cancer"]
    #             data.setdefault((
    #                     name,
    #                     treatment,
    #                     cell_line,
    #                     cancer
    #                 ), []
    #             )
    #             data[(
    #                     name,
    #                     treatment,
    #                     cell_line,
    #                     cancer
    #                 )
    #             ].append(float(line[i]))
    #         # For each sample...
    #         for s in data:
    #             # Initialize
    #             name = s[0]
    #             treatment = s[1]
    #             cell_line = s[2]
    #             cancer = s[3]
    #             # Get sample
    #             sample = Sample()
    #             if sample.is_unique(
    #                 session,
    #                 name,
    #                 treatment,
    #                 cell_line,
    #                 cancer
    #             ):
    #                 sample.name = name
    #                 sample.treatment = treatment
    #                 sample.cell_line = cell_line
    #                 sample.cancer = cancer
    #                 session.add(sample)
    #                 session.commit()
    #             sam = sample.select_unique(
    #                 session,
    #                 name,
    #                 treatment,
    #                 cell_line,
    #                 cancer
    #             )
    #             # Skip if feature not expressed in sample
    #             avg_expression_level = float(
    #                 sum(data[
    #                     name,
    #                     treatment,
    #                     cell_line,
    #                     cancer
    #                 ]) / len(data[
    #                     name,
    #                     treatment,
    #                     cell_line,
    #                     cancer
    #                 ])
    #             )
    #             if avg_expression_level > 0:
    #                 sampleIDs.append(sam.uid)
    #                 avg_expression_levels.append(
    #                     "%.3f" % avg_expression_level
    #                 )
    #         # Get TSS
    #         if feat_type == "tss":
    #             # Initialize
    #             tss_ids = set()
    #             for gene, tss_id in peak_ids:
    #                 tss = TSS()
    #                 if tss.is_unique(
    #                     session,
    #                     reg.uid,
    #                     sou.uid,
    #                     exp.uid,
    #                     gene,
    #                     tss_id
    #                 ):
    #                     tss.regionID = reg.uid
    #                     tss.sourceID = sou.uid
    #                     tss.experimentID = exp.uid
    #                     tss.gene = gene
    #                     tss.tss = tss_id
    #                     tss.sampleIDs = "{},".format(
    #                         ",".join(map(str, sampleIDs))
    #                     )
    #                     tss.avg_expression_levels = "{},".format(
    #                         ",".join(avg_expression_levels)
    #                     )
    #                     session.add(tss)
    #                     session.commit()
    #                 tss = tss.select_unique(
    #                     session,
    #                     reg.uid,
    #                     sou.uid,
    #                     exp.uid,
    #                     gene,
    #                     tss_id
    #                 )
    #                 tss_ids.add(tss)
    #         # For each sample...
    #         for i in range(len(sampleIDs)):
    #             if feat_type == "enhancer":
    #                 enhancer = Enhancer()
    #                 if enhancer.is_unique(
    #                     session,
    #                     reg.uid,
    #                     sou.uid,
    #                     sampleIDs[i],
    #                     exp.uid
    #                 ):
    #                     enhancer.regionID = reg.uid
    #                     enhancer.sourceID = sou.uid
    #                     enhancer.sampleID = sampleIDs[i]
    #                     enhancer.experimentID = exp.uid
    #                     # Append to features
    #                     features.append(enhancer)
    #             if feat_type == "tss":
    #                 for tss in tss_ids:
    #                     expression = Expression()
    #                     if expression.is_unique(
    #                         session,
    #                         tss.uid,
    #                         sampleIDs[i]
    #                     ):
    #                         expression.tssID = tss.uid
    #                         expression.sampleID = sampleIDs[i]
    #                         expression.avg_expression_level =\
    #                             avg_expression_levels[i]
    #                         # Append to features
    #                         features.append(expression)
    #         # Insert features
    #         session.add_all(features)
    #         session.commit()

#-------------#
# Main        #
#-------------#

if __name__ == "__main__": main()