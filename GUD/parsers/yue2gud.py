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

    def __init__(self, file_name, sample_name, experiment_type, restriction_enzyme, treatment, cell_line, cancer, sex, add):

        self.file_name = file_name
        self.sample_name = sample_name
        self.experiment_type = experiment_type
        self.restriction_enzyme = restriction_enzyme
        self.treatment = treatment
        self.cell_line = cell_line
        self.cancer = cancer
        self.sex = sex
        self.add = add

    @property
    def X(self):
        """
        number of X chromosomes
        """
        if self.sex == "female":
            return(2)
        elif self.sex == "male":
            return(1)
        else:
            return(None)

    @property
    def Y(self):
        """
        number of Y chromosomes
        """
        if self.sex == "female":
            return(0)
        elif self.sex == "male":
            return(1)
        else:
            return(None)

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
    e.g. python -m GUD.parsers.yue2gud --genome hg38 --samples ./samples/3dGenomeBrowser.tsv --feature tad --test -P 3306
    """

    # Initialize
    global chroms
    global engine
    global Feature
    global Session

    # Testing
    if test:
        global current_process
        from multiprocessing import current_process

    # Create dummy dir
    subdir = "%s.%s" % (genome, os.path.basename(__file__))
    dummy_dir = os.path.join(dummy_dir, subdir)
    if not os.path.isdir(dummy_dir):
        os.makedirs(dummy_dir)

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

    # Download data
    data_file, url = _download_data(genome, feat_type, dummy_dir)

    # Start a new session
    session = Session()

    # Get valid chromosomes
    chroms = ParseUtils.get_chroms(session)

    # Get source
    source = Source()
    source.name = "3D Genome Browser"
    source.url = url
    ParseUtils.upsert_source(session, source)
    source = ParseUtils.get_source(session, source.name, source.source_metadata, source.metadata_descriptor, url)

    # This is ABSOLUTELY necessary to prevent MySQL from crashing!
    session.close()
    engine.dispose()

    # Get samples
    samples = _get_samples(session, samples_file)

    # Extract data
    extracted_files = _extract_data(data_file, samples, dummy_dir)

    # For each file...
    for extracted_file in sorted(extracted_files):

        # Initialize
        extracted_file = extracted_file.split("/")
        tdgenbrow = samples[extracted_file[-1]]
        extracted_file = "/".join(extracted_file)

        # Split data
        data_files = _split_data(extracted_file, threads)

        # Skip        
        if tdgenbrow.add:

            # Initialize
            session = Session()

            # Get experiment
            experiment = Experiment()
            experiment.name = tdgenbrow.experiment_type
            experiment.experiment_metadata = "%s," % tdgenbrow.restriction_enzyme
            experiment.metadata_descriptor = "restriction_enzyme,"
            ParseUtils.upsert_experiment(session, experiment)
            experiment = ParseUtils.get_experiment(session, experiment.name, experiment.experiment_metadata, experiment.metadata_descriptor)

            # Get sample
            sample = Sample()
            sample.name = tdgenbrow.sample_name
            sample.treatment = tdgenbrow.treatment
            sample.cell_line = tdgenbrow.cell_line
            sample.cancer = tdgenbrow.cancer
            if tdgenbrow.sex is not None:
                sample.X = tdgenbrow.X
                sample.Y = tdgenbrow.Y
            ParseUtils.upsert_sample(session, sample)
            sample = ParseUtils.get_sample(session, sample.name, sample.X, sample.Y, sample.treatment, sample.cell_line, sample.cancer)

            # This is ABSOLUTELY necessary to prevent MySQL from crashing!
            session.close()
            engine.dispose()

            # Parallelize inserts to the database
            ParseUtils.insert_data_files_in_parallel(data_files, partial(_insert_data, experiment=experiment, sample=sample, source=source, test=test), threads)

    # Remove files
    if remove:
        shutil.rmtree(dummy_dir)

    # Remove session
    Session.remove()

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
        if line[7] != "male" and line[7] != "female":
            sex = None
        else:
            sex = line[7]
        if line[8] == "Yes":
            add = True
        else:
            add = False
        samples.setdefault(line[0], TDGenBrow(line[0], line[1], line[2], line[3], treatment, cell_line, cancer, sex, add))

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

def _insert_data(data_file, experiment, sample, source, test=False):

    # Initialize
    session = Session()

    # Testing
    if test:
        lines = 0
        print(current_process().name)

    # For each line...
    for line in ParseUtils.parse_tsv_file(data_file):

        # Skip empty lines
        if len(line) == 0:
            continue

        if Feature.__tablename__ == "tads":

            # Upsert region
            region = Region()
            region.chrom = str(line[0])
            if region.chrom.startswith("chr"):
                region.chrom = region.chrom[3:]
            if region.chrom not in chroms:
                continue
            region.start = int(line[1])
            region.end = int(line[2])
            region.bin = assign_bin(region.start, region.end)
            ParseUtils.upsert_region(session, region)

            # Get region
            region = ParseUtils.get_region(session, region.chrom, region.start, region.end)

            # Upsert feature
            feature = Feature()
            feature.region_id = region.uid
            feature.sample_id = sample.uid
            feature.experiment_id = experiment.uid
            feature.source_id = source.uid
            ParseUtils.upsert_tad(session, feature)

        else:
            pass

        # Testing
        if test:
            lines += 1
            if lines == 1000:
                break

    # This is ABSOLUTELY necessary to prevent MySQL from crashing!
    session.close()
    engine.dispose()

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":
    main()
