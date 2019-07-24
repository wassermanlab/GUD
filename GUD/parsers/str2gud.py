#!/usr/bin/env python

from . import ParseUtils
from GUD.ORM.source import Source
from GUD.ORM.region import Region
from GUD.ORM.short_tandem_repeat import ShortTandemRepeat
from GUD import GUDUtils
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

usage_msg = """
usage: %s --genome STR [-h] [options]
""" % os.path.basename(__file__)

help_msg = """%s
inserts short tandem repeats from curated lists to GUD.

  --genome STR      genome assembly
  --str_file FILE        short tandem repeat file   
  --based INT       1 or 0 based  
  --source STR      source name

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
    parser.add_argument("--str_file")
    parser.add_argument("--based")

    # Optional args
    optional_group = parser.add_argument_group("optional arguments")
    optional_group.add_argument("-h", "--help", action="store_true")
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
    if not args.genome or not args.str_file or not args.based or not args.source:
        error = ["%s\n%s" % (usage_msg, os.path.basename(
            __file__)), "error", "argument \"--genome\" \"--str\" \"--based\" \"--source\" is required\n"]
        print(": ".join(error))
        exit(0)

    # Check based is 0 or 1 
    try:
        args.based = int(args.based)
    except: 
        error = ["%s\n%s" % (usage_msg, os.path.basename(__file__)), "error",
                 "argument \"-t\" \"--based\"", "invalid int value", "\"%s\"\n" % args.based]
        print(": ".join(error))
        exit(0)

    if int(args.based) not in [0,1]:
        error = ["%s\n%s" % (usage_msg, os.path.basename(
            __file__)), "error", "argument \"--based\" must be 0 or 1\n"]
        print(": ".join(error))
        exit(0)

    # Check file exists
    if os.path.isfile(args.str_file) :
        error = ["%s\n%s" % (usage_msg, os.path.basename(
            __file__)), "error", "argument \"--str\" must be a valid file\n"]
        print(": ".join(error))
        exit(0)

    # Check "-t" argument
    try:
        args.threads = int(args.threads)
    except:
        error = ["%s\n%s" % (usage_msg, os.path.basename(__file__)), "error",
                 "argument \"-t\" \"--threads\"", "invalid int value", "\"%s\"\n" % args.threads]
        print(": ".join(error))
        exit(0)

    # Check MySQL password
    if not args.pwd:
        args.pwd = ""

    # Check MySQL port
    try:
        args.port = int(args.port)
    except:
        error = ["%s\n%s" % (usage_msg, os.path.basename(__file__)), "error",
                 "argument \"-P\" \"--port\"", "invalid int value", "\"%s\"\n" % args.port]
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
    str_to_gud(args.genome, args.str_file, args.based, args.source, args.threads)


def str_to_gud(genome, str_file, based, sourceName, threads=1):
    """
    python -m GUD.parsers.strgene2gud --genome hg19 --str_file <file> --based 0 --source <source>
    """

    # Globals
    global chroms
    global engine
    global Session
    global source

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
    ParseUtils.create_table(ShortTandemRepeat)

    # Start a new session
    session = Session()

    # Get valid chromosomes
    chroms = ParseUtils.get_chroms(session)

    # Get source
    source = Source()
    source.name = sourceName
    ParseUtils.upsert_source(session, source)
    source = ParseUtils.get_source(session, sourceName)

    # This is ABSOLUTELY necessary to prevent MySQL from crashing!
    session.close()
    engine.dispose()

    # Parallelize inserts to the database
    ParseUtils.process_data_in_chunks(
        str_file, _insert_data_in_chunks, threads)

    # Remove session
    Session.remove()


## pick up here 
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
        region = ParseUtils.get_region(
            session, region.chrom, region.start, region.end, region.strand)

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


if __name__ == "__main__":
    main()
