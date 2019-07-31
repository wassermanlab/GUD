#!/usr/bin/env python

from . import ParseUtils
from GUD.ORM.clinvar import ClinVar
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

usage_msg = """
usage: %s --genome STR [-h] [options]
""" % os.path.basename(__file__)

help_msg = """%s
inserts copy number variants from curated sourcesq.

  --genome STR        genome assembly
  --clinvar_file FILE Clinvar file with columns #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tN
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
    parser.add_argument("--clinvar_file")
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
    if not args.genome or not args.clinvar_file or not args.source_name:
        error = ["%s\n%s" % (usage_msg, os.path.basename(__file__)), "error",
                 "argument \"--genome\" \"--clinvar_file\" \"--source_name\"  are required\n"]
        print(": ".join(error))
        exit(0)

    # Check file exists
    if not os.path.isfile(args.clinvar_file):
        error = ["%s\n%s" % (usage_msg, os.path.basename(
            __file__)), "error", "argument \"--clinvar_file\" must be a valid existing file\n"]
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
    clinvar_to_gud(args.genome, args.source_name, args.clinvar_file,
                  args.test, args.threads)


def clinvar_to_gud(genome, source_name, clinvar_file, test=False, threads=1):
    """
    python -m GUD.parsers.str2gud --genome hg19 --source_name <name> --snv_file <FILE> 
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
    engine, Session = GUDUtils._get_engine_session(db_name)

    # Initialize parser utilities
    ParseUtils.genome = genome
    ParseUtils.dbname = db_name
    ParseUtils.engine = engine

    # If database does not exist...
    ParseUtils.initialize_gud_db()

    # Create table
    ParseUtils.create_table(ClinVar)

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
    data_file = clinvar_file

    # Split data
    data_files = _split_data(data_file, threads)

    # Parallelize inserts to the database
    ParseUtils.insert_data_files_in_parallel(
        data_files, partial(_insert_data, test=test), threads)

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
            cmd = 'parallel -j %s --pipe --block 2M -k grep "^%s[[:space:]]" < %s > %s' % (
                threads, chrom, data_file, split_file)
            subprocess.call(cmd, shell=True)

        # Append split file
        statinfo = os.stat(split_file)
        if statinfo.st_size:
            split_files.append(split_file)
        else:
            os.remove(split_file)

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

        # Get region
        region = Region()
        region.chrom = line[0]
        region.start = int(line[1]) -1
        region.end = region.start + len(line[3])
        region.bin = assign_bin(region.start, region.end)

        # Ignore non-standard chroms, scaffolds, etc.
        if region.chrom not in chroms:
            continue

        # Upsert region
        ParseUtils.upsert_region(session, region)

        # Get region ID
        region = ParseUtils.get_region(
            session, region.chrom, region.start, region.end, region.strand)

        # Get feature
        clinvar = ClinVar()
        clinvar.region_id = region.uid
        clinvar.source_id = source.uid
        clinvar.ref = line[3].encode(encoding='UTF-8')
        clinvar.alt = line[4].encode(encoding='UTF-8')
        clinvar.clinvar_variation_ID = line[2]
        ## get info 
        info = line[7].split(";")
        infoDict = {}
        for t in info: 
            if "=" in t:
                infoDict[t.split("=")[0]] = t.split("=")[1]
        try:
            ANN = infoDict["ANN"].split("|")
            clinvar.ANN_Annotation = ANN[1].encode(encoding='UTF-8')
            clinvar.ANN_Annotation_Impact = ANN[2].encode(encoding='UTF-8')
            clinvar.ANN_Gene_Name = ANN[3].encode(encoding='UTF-8')
            clinvar.ANN_Gene_ID = ANN[4].encode(encoding='UTF-8')
            clinvar.ANN_Feature_Type = ANN[5].encode(encoding='UTF-8')
            clinvar.ANN_Feature_ID = ANN[6].encode(encoding='UTF-8')
        except:
            clinvar.ANN_Annotation = None
            clinvar.ANN_Annotation_Impact = None
            clinvar.ANN_Gene_Name = None
            clinvar.ANN_Gene_ID = None
            clinvar.ANN_Feature_Type = None
            clinvar.ANN_Feature_ID = None
        try: 
            clinvar.CADD = float(infoDict["CADD"])
        except:  
            clinvar.CADD = None
        try:
            clinvar.CLNDISDB = infoDict["CLNDISDB"].encode(encoding='UTF-8')
        except:  
            clinvar.CLNDN = None
        try:
            clinvar.CLNDN = infoDict["CLNDN"].encode(encoding='UTF-8')
        except:  
            clinvar.CLNDN = None
        try:
            clinvar.CLNSIG = infoDict["CLNSIG"].encode(encoding='UTF-8')
        except:  
            clinvar.CLNSIG = None
        try:
            clinvar.gnomad_exome_af_global      = float(infoDict["gnomad_exome_af_global"])
        except:
            clinvar.gnomad_exome_af_global      = None
        try:
            clinvar.gnomad_exome_hom_global     = float(infoDict["gnomad_exome_hom_global"])
        except:
            clinvar.gnomad_exome_hom_global     = None
        try:
            clinvar.gnomad_genome_af_global     = float(infoDict["gnomad_genome_af_global"])
        except:
            clinvar.gnomad_genome_af_global     = None
        try:
            clinvar.gnomad_genome_hom_global    = float(infoDict["gnomad_genome_hom_global"])
        except:
            clinvar.gnomad_genome_hom_global    = None

        ParseUtils.upsert_clinvar(session, clinvar)
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