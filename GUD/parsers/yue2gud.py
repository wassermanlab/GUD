#!/usr/bin/env python

import argparse
from binning import assign_bin
from functools import partial
import getpass
from multiprocessing import Pool, cpu_count
from numpy import isnan
import os
import pickle
import re
import requests
import shutil
import subprocess
import sys
import warnings
from zipfile import ZipFile
# Python 3+
if sys.version_info > (3, 0):
    from urllib.request import urlretrieve
# Python 2.7
else:
    from urllib import urlretrieve

# Import from GUD module
from GUD import GUDUtils
from GUD.ORM.experiment import Experiment
from GUD.ORM.region import Region
from GUD.ORM.sample import Sample
from GUD.ORM.source import Source
from GUD.ORM.tad import TAD
from . import ParseUtils

usage_msg = """
usage: %s --genome STR --samples FILE --feature STR [-h] [options]
""" % os.path.basename(__file__)

help_msg = """%s
inserts features from the 3D Genome Browser, which is developed by
the Yue lab, into GUD.

  --genome STR        genome assembly
  --samples FILE      3D Genome Browser samples (manually-curated)
  --feature STR       type of genomic feature ("loop" or "tad")

optional arguments:
  -h, --help          show this help message and exit
  --dummy-dir DIR     dummy directory (default = "/tmp/")
  -r, --remove        remove downloaded files (default = False)
  -t, --test          limit the total of inserts to ~1K per
                      thread for testing (default = False)
  --threads INT       number of threads to use (default = %s)

mysql arguments:
  -d STR, --db STR    database name (default = "%s")
  -H STR, --host STR  host name (default = "localhost")
  -p STR, --pwd STR   password (default = ignore this option)
  -P INT, --port INT  port number (default = %s)
  -u STR, --user STR  user name (default = current user)
""" % (usage_msg, (cpu_count() - 1), GUDUtils.db, GUDUtils.port)

#-------------#
# Classes     #
#-------------#

class TDGenBrow:

    def __init__(self, file_name, sample_name, experiment_type, restriction_enzyme, treatment, cell_line, cancer, add):

        self.file_name = file_name
        self.sample_name = sample_name
        self.experiment_type = experiment_type
        self.restriction_enzyme = restriction_enzyme
        self.treatment = treatment
        self.cell_line = cell_line
        self.cancer = cancer
        self.add = add

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
    parser.add_argument("--feature")

    # Optional args
    optional_group = parser.add_argument_group("optional arguments")
    optional_group.add_argument("-h", "--help", action="store_true")
    optional_group.add_argument("--dummy-dir", default="/tmp/")
    optional_group.add_argument("-r", "--remove", action="store_true")
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

    # Check for invalid feature
    valid_features = ["loop", "tad"]
    if args.feature not in valid_features:
        error = ["%s\n%s" % (usage_msg, os.path.basename(__file__)), "error", "argument \"--feature\"", "invalid choice", "\"%s\" (choose from" % args.feature, "%s)\n" % " ".join(["\"%s\"" % i for i in valid_features])]
        print(": ".join(error))
        exit(0)

    # Check "--threads" argument
    try:
        args.threads = int(args.threads)
    except:
        error = ["%s\n%s" % (usage_msg, os.path.basename(__file__)), "error", "argument \"-t\" \"--threads\"", "invalid int value", "\"%s\"\n" % args.threads]
        print(": ".join(error))
        exit(0)

    # Check MySQL password
    if not args.pwd:
        args.pwd = ""

    # Check MySQL port
    try:
        args.port = int(args.port)
    except:
        error = ["%s\n%s" % (usage_msg, os.path.basename(__file__)), "error", "argument \"-P\" \"--port\"", "invalid int value", "\"%s\"\n" % args.port]
        print(": ".join(error))
        exit(0)

def main():

    # Parse arguments
    args = parse_args()

    # Set MySQL options
    GUDUtils.user = args.user
    GUDUtils.pwd = args.pwd
    GUDUtils.host = args.host
    GUDUtils.port = args.port
    GUDUtils.db = args.db

    # Insert data from the 3D Genome Browser
    yue_to_gud(args.genome, args.samples, args.feature, args.dummy_dir, args.remove, args.test, args.threads)

def yue_to_gud(genome, samples_file, feat_type, dummy_dir="/tmp/", remove=False, test=False, threads=1):
    """
    e.g. python -m GUD.parsers.yue2gud --genome hg38 --feature abcd --test -P 3306
    """

    # Initialize
    global chroms
    global engine
    global experiment
    global restriction_enzyme
    global sample
    global source
    global Feature
    global Session

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
    if feat_type == "tad":
        Feature = TAD
    else:
        pass
    ParseUtils.create_table(Feature)

    # Start a new session
    session = Session()

    # Get valid chromosomes
    chroms = ParseUtils.get_chroms(session)

    # Download data
    data_file, url = _download_data(genome, feat_type, dummy_dir)

    # Get source
    source = Source()
    source.name = "3D Genome Browser"
    source.url = url
    ParseUtils.upsert_source(session, source)
    source = ParseUtils.get_source(session, source.name, source.source_metadata, source.metadata_descriptor, url)

    # Get samples
    samples = _get_samples(session, samples_file)

    # This is ABSOLUTELY necessary to prevent MySQL from crashing!
    session.close()
    engine.dispose()

    # Extract data
    extracted_files = _extract_data(data_file, samples, dummy_dir)

    # For each file...
    for extracted_file in sorted(extracted_files):

        # Split data
        data_files = _split_data(extracted_file, threads)

        # Skip
        extracted_file = extracted_file.split("/")
        if not samples[extracted_file[-1]].add:
            continue

        # Get experiment
        experiment = Experiment()
        experiment.name = samples[extracted_file[-1]].experiment_type
        experiment.experiment_metadata = "%s," % samples[extracted_file[-1]].restriction_enzyme
        experiment.metadata_descriptor = "restriction_enzyme,"
        ParseUtils.upsert_experiment(session, experiment)
        experiment = ParseUtils.get_experiment(session, experiment.name, experiment.experiment_metadata, experiment.metadata_descriptor)

        # Get sample
        sample = Sample()
        sample.name = samples[extracted_file[-1]].sample_name
        sample.treatment = samples[extracted_file[-1]].sample_name
        sample.cell_line = samples[extracted_file[-1]].sample_name
        sample.cancer = samples[extracted_file[-1]].sample_name
        if encodes[accession].sex is not None:
            sample.X = encodes[accession].X
            sample.Y = encodes[accession].Y       
        ParseUtils.upsert_sample(session, sample)


        self.file_name = file_name
        self.sample_name = sample_name
        self.experiment_type = experiment_type
        self.restriction_enzyme = restriction_enzyme
        self.treatment = treatment
        self.cell_line = cell_line
        self.cancer = cancer
        self.add = add

        # # Parallelize inserts to the database
        # ParseUtils.insert_data_files_in_parallel(data_files, partial(_insert_data, test=test), threads)

        # # Remove data files
        # if remove:
        #     if os.path.exists(data_file):
        #         os.remove(data_file)
        #     for data_file in data_files:
        #         if os.path.exists(data_file):
        #             os.remove(data_file)

        # # Remove session
        # Session.remove()

def _download_data(genome, feat_type, dummy_dir="/tmp/"):

    # Initialize
    url = "http://promoter.bx.psu.edu/hi-c/downloads/"

    # Get ftp file
    if feat_type == "loops":
        ftp_file = "loops-%s.zip" % genome
    else:
        ftp_file = "%s.TADs.zip" % genome

    # Download data
    data_file = os.path.join(dummy_dir, ftp_file)
    if not os.path.exists(data_file):
        f = urlretrieve(os.path.join(url, ftp_file), data_file)

    return(data_file, os.path.join(url, ftp_file))

def _get_samples(session, file_name):

    # Initialize
    samples = {}

    # For each line...
    for line in ParseUtils.parse_tsv_file(file_name):

        # Add 3D Genome Browser object
        if line[4] == "Yes":
            treatment = True
        else:
            treatment = False
        if line[5] == "Yes":
            cell_line = True
        else:
            cell_line = False
        if line[6] == "Yes":
            cancer = True
        else:
            cancer = False
        if line[7] == "Yes":
            add = True
        else:
            add = False
        samples.setdefault(line[0], TDGenBrow(line[0], line[1], line[2], line[3], treatment, cell_line, cancer, add))

    return(samples)

def _extract_data(data_file, samples, dummy_dir):

    # Initialize
    extracted_files = []
    i_have_been_warned = False

    # For each data file...
    with ZipFile(data_file) as myzip:

        for zip_file in myzip.namelist():

            zip_file = zip_file.split("/")

            # Skip
            if zip_file[-1] == "":
                continue

            # Warn me!
            if zip_file[-1] not in samples:
                warnings.warn("missing sample: %s" % zip_file[-1], Warning, stacklevel=2)
                i_have_been_warned = True

            # Extract file
            extracted_file = os.path.join(dummy_dir, zip_file[-1])
            if not os.path.exists(extracted_file):
                content = myzip.read("/".join(zip_file))
                ParseUtils.write(extracted_file, content.decode(encoding="UTF-8"))

            # Add extracted file
            extracted_files.append(extracted_file)

    if i_have_been_warned:
        error = ["%s" % os.path.basename(__file__), "error", "missing samples"]
        print(": ".join(error))
        exit(0)

    return(extracted_files)

def _split_data(data_file, threads=1):

    # Initialize
    split_files = []
    split_dir = os.path.dirname(os.path.realpath(data_file))

    # Get number of lines
    output = subprocess.check_output(["zless %s | wc -l" % data_file], shell=True)
    m = re.search("(\d+)", str(output))
    L = float(m.group(1))

    # Split
    prefix = "%s." % data_file.split("/")[-1]
    cmd = "zless %s | split -d -l %s - %s" % (data_file, int(L/threads)+1, os.path.join(split_dir, prefix))
    subprocess.run(cmd, shell=True)

    # For each split file...
    for split_file in os.listdir(split_dir):

        # Append split file
        if split_file.startswith(prefix):
            split_files.append(os.path.join(split_dir, split_file))

    return(split_files)

def _insert_data(data_file, test=False):

    # Initialize
    session = Session()

    # Testing
    if test:
        lines = 0
        print(current_process().name)

    # For each line...
    for line in ParseUtils.parse_tsv_file(data_file):

        # Skip empty lines
        if not line:
            continue

        # Upsert region
        region = Region()
        region.chrom = line[5][3:]
        region.start = line[6]
        region.end = line[7]
        region.bin = assign_bin(region.start, region.end)
        if region.chrom not in chroms:
            continue
        ParseUtils.upsert_region(session, region)

        # Get region
        region = ParseUtils.get_region(session, region.chrom, region.start, region.end)

        # Upsert CpG island
        repeat = RepeatMask()
        repeat.score = line[1]
        repeat.name = line[10]
        repeat.repeat_class = line[11]
        repeat.family = line[12]
        repeat.strand = line[9]
        repeat.region_id = region.uid
        repeat.source_id = source.uid
        ParseUtils.upsert_rmsk(session, repeat)

        # Testing
        if test:
            lines += 1
            if lines > 1000:
                break

    # This is ABSOLUTELY necessary to prevent MySQL from crashing!
    session.close()
    engine.dispose()

#-------------#
# Main        #
#-------------#

if __name__ == "__main__": main()