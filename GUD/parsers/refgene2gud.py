#!/usr/bin/env python

import argparse
from binning import assign_bin
import getpass
from multiprocessing import cpu_count
from multiprocessing import current_process
import os
import re
import sys
# Python 3+
if sys.version_info > (3, 0):
    from urllib.request import urlretrieve
# Python 2.7
else:
    from urllib import urlretrieve

# Import from GUD module
from GUD import GUDUtils
from GUD.ORM.gene import Gene
from GUD.ORM.region import Region
from GUD.ORM.source import Source
from . import ParseUtils

usage_msg = """
usage: %s --genome STR [-h] [options]
""" % os.path.basename(__file__)

help_msg = """%s
inserts features from the UCSC's "refGene" table into GUD.

  --genome STR        genome assembly

optional arguments:
  -h, --help          show this help message and exit
  --dummy-dir DIR     dummy directory (default = "/tmp/")
  --threads INT       number of additional threads to use
                      (default = %s)

mysql arguments:
  -d STR, --db STR    database name (default = "%s")
  -H STR, --host STR  host name (default = "localhost")
  -p STR, --pwd STR   password (default = ignore this option)
  -P INT, --port INT  port number (default = %s)
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

    # Optional args
    optional_group = parser.add_argument_group("optional arguments")
    optional_group.add_argument("-h", "--help", action="store_true")
    optional_group.add_argument("--dummy-dir", default="/tmp/")
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
    if not args.genome:
        error = ["%s\n%s" % (usage_msg, os.path.basename(__file__)), "error", "argument \"--genome\" is required\n"]
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

    # Insert RefGene data
    refgene_to_gud(args.genome, args.dummy_dir, args.threads)

def refgene_to_gud(genome, dummy_dir="/tmp/", threads=1):
    """
    python -m GUD.parsers.refgene2gud --genome hg19 --dummy-dir ./tmp/
    """

    # Globals
    global chroms
    global engine
    global Session
    global source

    # Download data
    dummy_file = _download_data(genome, dummy_dir)

    # Get database name
    db_name = GUDUtils._get_db_name()

    # Get engine/session
    engine, Session = GUDUtils._get_engine_session(db_name)

    # Initialize parser utilities
    ParseUtils.genome = genome
    ParseUtils.dbname = db_name
    ParseUtils.engine = engine

    # If database does not exist...
    ParseUtils.initialize_gud_db()

    # Create table
    ParseUtils.create_table(Gene)

    # Start a new session
    session = Session()

    # Get valid chromosomes
    chroms = ParseUtils.get_chroms(session)

    # Get source
    source = Source()
    m = re.search("^%s/*(.+).txt.gz$" % dummy_dir, dummy_file) 
    source_name = m.group(1)
    source.name = source_name
    ParseUtils.upsert_source(session, source)
    source = ParseUtils.get_source(session, source_name)

    # This is ABSOLUTELY necessary to prevent MySQL from crashing!
    session.close()
    engine.dispose()

    # Parallelize inserts to the database
    ParseUtils.process_data_in_chunks(dummy_file, _insert_data_in_chunks, threads)

    # Remove session
    Session.remove()

    # # Remove downloaded file
    # if os.path.exists(dummy_file):
    #     os.remove(dummy_file)

def _download_data(genome, dummy_dir="/tmp/"):

    # Initialize
    url = "ftp://hgdownload.soe.ucsc.edu/goldenPath/%s/database" % genome
    ftp_file = os.path.join(url, "refGene.txt.gz")
    dummy_file = os.path.join(dummy_dir, "refGene.txt.gz")

    # Download data
    if not os.path.exists(dummy_file):
        f = urlretrieve(ftp_file, dummy_file)

    return(dummy_file)

def _insert_data_in_chunks(chunk):

    print(current_process().name)

    # Start a new session
    session = Session()

    # For each line...
    for line in chunk:

        # Skip empty lines
        if not line:
            continue

        # Get region
        region = Region()
        region.chrom = line[2]
        region.start = int(line[4])
        region.end = int(line[5])
        region.bin = assign_bin(region.start, region.end)
        region.strand = line[3]

        # Ignore non-standard chroms, scaffolds, etc.
        if region.chrom not in chroms:
            continue

        # Upsert region
        ParseUtils.upsert_region(session, region)

        # Get region ID
        region = ParseUtils.get_region(session, region.chrom, region.start, region.end, region.strand)

        # Get gene
        gene = Gene()
        gene.regionID = region.uid
        gene.name = line[1]
        gene.name2 = line[12]
        gene.cdsStart = int(line[6])
        gene.cdsEnd = int(line[7])
        gene.exonStarts = line[9].encode("utf-8")
        gene.exonEnds = line[10].encode("utf-8")
        gene.sourceID = source.uid

        # Upsert gene
        ParseUtils.upsert_gene(session, gene)

    # This is ABSOLUTELY necessary to prevent MySQL from crashing!
    session.close()
    engine.dispose()

#-------------#
# Main        #
#-------------#

if __name__ == "__main__": main()