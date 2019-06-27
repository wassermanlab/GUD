#!/usr/bin/env python

import argparse
from binning import assign_bin
import getpass
from multiprocessing import cpu_count
from multiprocessing import current_process
import os
import re
from sqlalchemy_utils import database_exists
import sys
# Python 3+
if sys.version_info > (3, 0):
    from urllib.request import urlretrieve
# Python 2.7
else:
    from urllib import urlretrieve

# Import from GUD module
from GUD import GUDglobals
from GUD.ORM.gene import Gene
from GUD.ORM.region import Region
from GUD.ORM.source import Source
from . import _get_chroms, _get_db_name, _get_region, _get_source, _initialize_gud_db, _initialize_engine_session, _process_data_in_chunks, _upsert_gene, _upsert_region, _upsert_source

usage_msg = """
usage: %s --genome STR [-h] [options]
""" % os.path.basename(__file__)

help_msg = """%s
inserts features from the UCSC's "refGene" table into GUD.

  --genome STR        genome assembly

optional arguments:
  -h, --help          show this help message and exit
  --dummy-dir DIR     dummy directory (default = "/tmp/")
  -t, --threads       number of additional threads to use
                      (default = %s)

mysql arguments:
  -d STR, --db STR    database name (default = "%s")
  -H STR, --host STR  host name (default = "localhost")
  -p STR, --pwd STR   password (default = ignore this option)
  -P STR, --port STR  port number (default = %s)
  -u STR, --user STR  user name (default = current user)
""" % (usage_msg, (cpu_count() - 1), GUDglobals.db_name, GUDglobals.db_port)

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
    optional_group.add_argument("-t", "--threads", default=(cpu_count() - 1))
    
    # MySQL args
    mysql_group = parser.add_argument_group("mysql arguments")
    mysql_group.add_argument("-d", "--db", default=GUDglobals.db_name)
    mysql_group.add_argument("-H", "--host", default="localhost")
    mysql_group.add_argument("-p", "--pwd")
    mysql_group.add_argument("-P", "--port", default=GUDglobals.db_port)
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

def main():

    # Parse arguments
    args = parse_args()

    # Insert RefGene data
    refgene_to_gud(args.user, args.pwd, args.host, args.port, args.db, args.genome, args.dummy_dir, args.threads)

def refgene_to_gud(user, pwd, host, port, db, genome, dummy_dir="/tmp/", threads=1):

    # Globals
    global chroms
    global engine
    global Session
    global source

    # Get database name
    db_name = _get_db_name(user, pwd, host, port, db)

    # If database does not exist...
    if not database_exists(db_name):
        _initialize_gud_db(user, pwd, host, port, db, genome)

    # Get engine/session
    engine, Session = _initialize_engine_session(db_name)
    session = Session()

    # Create table
    if not engine.has_table(Gene.__tablename__):
        Gene.__table__.create(bind=engine)

    # Get valid chromosomes
    chroms = _get_chroms(session)

    # Download data
    dummy_file = _download_data(genome, dummy_dir)

    # Get source
    source = Source()
    m = re.search("^%s/*(.+).txt.gz$" % dummy_dir, dummy_file) 
    source_name = m.group(1)
    source.name = source_name
    _upsert_source(session, source)
    source = _get_source(session, source_name)

    # Close session
    session.close()

    # Dispose of session
    engine.dispose()

    # Parallelize inserts to the database
    _process_data_in_chunks(dummy_file, _insert_data_in_chunks, threads)

    # Dispose session
    Session.remove()

    # Remove downloaded file
    if os.path.exists(dummy_file):
        os.remove(dummy_file)

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

    # Initialize
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
        _upsert_region(session, region)

        # Get region ID
        region = _get_region(session, region.chrom, region.start, region.end, region.strand)

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
        _upsert_gene(session, gene)

    # Close session
    session.close()

    # Dispose of session
    engine.dispose()

#-------------#
# Main        #
#-------------#

if __name__ == "__main__": main()