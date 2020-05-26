#!/usr/bin/env python

import argparse
from binning import assign_bin
import getpass
import os
import re
import shutil
import subprocess
import sys
import warnings

# Import from GUD module
from GUD import GUDUtils
from GUD.ORM.enhancer import Enhancer
from GUD.ORM.experiment import Experiment
from GUD.ORM.region import Region
from GUD.ORM.sample import Sample
from GUD.ORM.source import Source
from . import ParseUtils

usage_msg = """
usage: %s --genome STR [-h] [options]
""" % os.path.basename(__file__)

help_msg = """%s
inserts positive enhancers from the VISTA into GUD.

  --genome STR        genome assembly

optional arguments:
  -h, --help          show this help message and exit
  --dummy-dir DIR     dummy directory (default = "/tmp/")
  -r, --remove        remove downloaded files (default = False)

mysql arguments:
  -d STR, --db STR    database name (default = "%s")
  -H STR, --host STR  host name (default = "localhost")
  -p STR, --pwd STR   password (default = ignore this option)
  -P INT, --port INT  port number (default = %s)
  -u STR, --user STR  user name (default = current user)
""" % (usage_msg, GUDUtils.db, GUDUtils.port)

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

    # Optional args
    optional_group = parser.add_argument_group("optional arguments")
    optional_group.add_argument("-h", "--help", action="store_true")
    optional_group.add_argument("--dummy-dir", default="/tmp/")
    optional_group.add_argument("-r", "--remove", action="store_true")
    
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
    if not args.genome:
        error = ["%s\n%s" % (usage_msg, os.path.basename(__file__)), "error", "argument \"--genome\" is required\n"]
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

    # Insert VISTA data
    vista_to_gud(args.genome, args.dummy_dir, args.remove)

def vista_to_gud(genome, dummy_dir="/tmp/", remove=False):
    """
    e.g. python -m GUD.parsers.vista2gud --genome hg38 -P 3306
    """

    # Initialize
    global chroms
    global engine
    global experiment
    global session
    global source

    # Create dummy dir
    subdir = "%s.%s" % (genome, os.path.basename(__file__))
    dummy_dir = os.path.join(dummy_dir, subdir)
    if not os.path.isdir(dummy_dir):
        os.makedirs(dummy_dir)

    # Download data
    data_files, url = _download_data(genome, dummy_dir)

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
    ParseUtils.create_table(Enhancer)

    # Start a new session
    session = Session()

    # Get valid chromosomes
    chroms = ParseUtils.get_chroms(session)

    # Get experiment
    experiment = Experiment()
    experiment.name = "transgenic reporter assay"
    ParseUtils.upsert_experiment(session, experiment)
    experiment = ParseUtils.get_experiment(session, experiment.name)    

    # Get source
    source = Source()
    source.name = "VISTA"
    if genome == "hg38" or genome == "mm10":
        source.source_metadata = "True,"
    else:
        source.source_metadata = "False,"
    source.metadata_descriptor = "liftOver,"
    source.url = url
    ParseUtils.upsert_source(session, source)
    source = ParseUtils.get_source(session, source.name, source.source_metadata,
                                   source.metadata_descriptor, url)

    # Insert to the database
    _insert_data(genome, data_files[0], data_files[1])

    # Remove files
    if remove:
        shutil.rmtree(dummy_dir)

    # This is ABSOLUTELY necessary to prevent MySQL from crashing!
    session.close()
    engine.dispose()

    # Remove session
    Session.remove()

def _download_data(genome, dummy_dir="/tmp/"):

    # Initialize
    dummy_files = []

    # Python 3+
    if sys.version_info > (3, 0):
        from urllib.request import urlretrieve
    # Python 2.7
    else:
        from urllib import urlretrieve

    if genome == "hg38" or genome == "mm10":

        if genome == "hg38":
            url = "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/"
            chains_file = "hg19ToHg38.over.chain.gz"
        else:
            url = "http://hgdownload.soe.ucsc.edu/goldenPath/mm9/liftOver/"
            chains_file = "mm9ToMm10.over.chain.gz"

        dummy_files.append(os.path.join(dummy_dir, chains_file))

        # Download data
        if not os.path.exists(dummy_files[-1]):
            f = urlretrieve(os.path.join(url, chains_file), dummy_files[-1])

    else:
        dummy_files.append(None)

    # Download data
    url = "https://enhancer.lbl.gov/cgi-bin/imagedb3.pl?page_size=20000;show=1;search.result=yes;page=1;form=search;search.form=no;action=search;search.sequence=1"
    dummy_files.append(os.path.join(dummy_dir, "vista.fa"))
    if not os.path.exists(dummy_files[-1]):
        f = urlretrieve(url, dummy_files[-1])

    return(dummy_files[::-1], url)

def _insert_data(genome, data_file, chains_file=None):

    # Initialize
    samples = {}
    if chains_file is not None:
        from pyliftover import LiftOver
        lo = LiftOver(chains_file)

    # For each SeqRecord...
    for seq_record in ParseUtils.parse_fasta_file(data_file, True):

        line = [s.strip() for s in seq_record.description.split("|")]

        if line[3] == "negative":
            continue

        if line[0] == "Human":
            if genome == "mm9" or genome == "mm10":
                continue
        else:
            if genome == "hg19" or genome == "hg38":
                continue

        # Upsert region
        region = Region()
        m = re.search("(chr\S+):(\d+)-(\d+)", line[1])
        region.chrom = m.group(1)
        if region.chrom.startswith("chr"):
            region.chrom = region.chrom[3:]
        if region.chrom not in chroms:
            continue
        region.start = int(m.group(2)) - 1 # i.e. 1-based
        region.end = int(m.group(3))
        if chains_file:
            try:
                chrom = "chr%s" % region.chrom
                start = lo.convert_coordinate(chrom, region.start)
                region.start = start[0][1]
                end = lo.convert_coordinate(chrom, region.end - 1) # i.e. requires 0-based
                region.end = end[0][1]
            except:
                msg = "position could not be found in new assembly"
                warnings.warn("%s: %s" % (msg, line[1]), Warning, stacklevel=2)
                continue
        region.bin = assign_bin(region.start, region.end)
        ParseUtils.upsert_region(session, region)

        # Get region
        region = ParseUtils.get_region(session, region.chrom, region.start, region.end)

        # For each sample...
        for s in line[4:]:

            if s not in samples:

                # Upsert sample
                sample = Sample()
                m = re.search("(.+)\[\d+\/\d+\]", s)
                sample.name = m.group(1)
                sample.treatment = False
                sample.cell_line = False
                sample.cancer = False
                ParseUtils.upsert_sample(session, sample)

                # Get sample
                sample = ParseUtils.get_sample(session, sample.name, sample.X, sample.Y, sample.treatment, sample.cell_line, sample.cancer)
                samples.setdefault(s, sample)
    
            # Upsert enhancer
            enhancer = Enhancer()
            enhancer.region_id = region.uid
            enhancer.experiment_id = experiment.uid
            enhancer.source_id = source.uid
            enhancer.sample_id = samples[s].uid
            ParseUtils.upsert_enhancer(session, enhancer)

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":
    main()
