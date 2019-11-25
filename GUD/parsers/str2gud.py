#!/usr/bin/env python

from . import ParseUtils
from GUD.ORM.short_tandem_repeat import ShortTandemRepeat
from GUD.ORM.source import Source
from GUD.ORM.region import Region
from GUD import GUDUtils
from binning import assign_bin
from functools import partial
import getpass
from multiprocessing import Pool, cpu_count
from numpy import isnan
import os
import sys
import re
import subprocess
import warnings
import argparse

# Import from GUD module
# python -m GUD.parsers.str2gud --test \
# --genome hg38 --source_name GangSTR \
# --str_file /space/home/tavshalom/GUD_TABLES/Transfer_GRCh38_190723/Genomic_STR.GUDformatted.tsv \
# --based 1 -d hg38 -u gud_w 
usage_msg = """
usage: %s --genome STR [-h] [options]
""" % os.path.basename(__file__)

help_msg = """%s
inserts short tandem repeats into GUD from curated sources
into GUD.

  --genome STR        genome assembly
  --str_file FILE     Short Tandem Repeat file with columns #chrom\tstart\tend\tmotif\tpathogenicity
  --based INT         0 or 1 based
  --source_name STR   source name      

optional arguments:
  -h, --help          show this help message and exit
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
    parser.add_argument("--source_name")

    # Optional args
    optional_group = parser.add_argument_group("optional arguments")
    optional_group.add_argument("-h", "--help", action="store_true")
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
    if not args.genome or not args.str_file or not args.based or not args.source_name:
        error = ["%s\n%s" % (usage_msg, os.path.basename(__file__)), "error",
                 "argument \"--genome\" \"--str_file\" \"--based\" \"--source_name\"  are required\n"]
        print(": ".join(error))
        exit(0)

    # Check file exists
    if not os.path.isfile(args.str_file):
        error = ["%s\n%s" % (usage_msg, os.path.basename(
            __file__)), "error", "argument \"--str_file\" must be a valid existing file\n"]
        print(": ".join(error))
        exit(0)

    # Check based
    try:
        args.based = int(args.based)
        if args.based not in [0, 1]:
            error = ["%s\n%s" % (usage_msg, os.path.basename(
                __file__)), "error", "argument \"--based\" must be 0 or 1\n"]
            print(": ".join(error))
            exit(0)
    except:
        error = ["%s\n%s" % (usage_msg, os.path.basename(
            __file__)), "error", "argument \"--based\" must be 0 or 1\n"]
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

# str_file, based, source_name


def main():

    # Parse arguments
    args = parse_args()

    # Set MySQL options
    GUDUtils.user = args.user
    GUDUtils.pwd = args.pwd
    GUDUtils.host = args.host
    GUDUtils.port = args.port
    GUDUtils.db = args.db

    # Insert ENCODE data
    str_to_gud(args.genome, args.source_name, args.str_file,
                  args.based, args.test, args.threads)


# TODO removed m and dummy_dir
def str_to_gud(genome, source_name, str_file, based, test=False, threads=1):
    """
    python -m GUD.parsers.str2gud --genome hg19 --source_name <name> --str_file <FILE> --based <0|1>
    """
    # Initialize
    global chroms
    global engine
    global Session
    global source

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
    ParseUtils.create_table(ShortTandemRepeat)

    # Start a new session
    session = Session()

    # Get valid chromosomes
    chroms = ParseUtils.get_chroms(session)


    # Get source
    source = Source()
    source.name = source_name
    ParseUtils.upsert_source(session, source)
    source = ParseUtils.get_source(session, source_name)

    # This is ABSOLUTELY necessary to prevent MySQL from crashing!
    session.close()
    engine.dispose()

    # Prepare data
    data_file = str_file

    # Split data
    data_files = _split_data(data_file, threads)

    # Parallelize inserts to the database
    ParseUtils.insert_data_files_in_parallel(
        data_files, partial(_insert_data, based=based, test=test), threads)

    # Remove data file
    for df in data_files:
        if os.path.exists(df):
            os.remove(df)

    # Remove session
    Session.remove()

def _split_data(data_file, threads=1):

    # import math

    # Initialize
    split_files = []

    # For each chromosome...
    for chrom in chroms:

        # Skip if file already split
        split_file = "%s.%s" % (data_file, chrom)
        if not os.path.exists(split_file):

            # Parallel split
            # cmd = 'parallel -j %s --pipe --block 2M -k grep "^%s[[:space:]]" < %s > %s' % (
            #     threads, chrom, data_file, split_file)
            cmd = 'parallel -j %s --pipe --block 2M -k grep "^%s[[:space:]]" < %s > %s' % (
                threads, "chr%s" % chrom, data_file, split_file)
            subprocess.call(cmd, shell=True)

        # Append split file
        statinfo = os.stat(split_file)
        if statinfo.st_size:
            split_files.append(split_file)
        else:
            os.remove(split_file)

    return(split_files)


def _insert_data(data_file, based=1, test=False):

    # Initialize
    session = Session()

    # Testing
    if test:
        lines = 0
        print(current_process().name)

    # For each line...
    for line in ParseUtils.parse_tsv_file(data_file):

        # Get region
        region = Region()
        # region.chrom = line[0]
        region.chrom = line[0][3:]
        region.start = int(line[1])
        if based:
            region.start -= 1
        region.end = int(line[2])
        region.bin = assign_bin(region.start, region.end)

        # Ignore non-standard chroms, scaffolds, etc.
        if region.chrom not in chroms:
            continue

        # Upsert region
        ParseUtils.upsert_region(session, region)

        # Get region ID
        region = ParseUtils.get_region(
            session, region.chrom, region.start, region.end)

        # Get feature
        STR = ShortTandemRepeat()
        STR.region_id = region.uid
        STR.source_id = source.uid
        STR.motif = line[3]
        STR.pathogenicity = int(line[4])

        ParseUtils.upsert_str(session, STR)

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