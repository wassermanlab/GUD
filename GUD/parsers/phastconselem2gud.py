#!/usr/bin/env python

import argparse
from binning import assign_bin
from functools import partial
import getpass
from multiprocessing import cpu_count
import os
import re
import shutil
import subprocess
import sys

# Import from GUD module
from GUD import GUDUtils
from GUD.ORM.conservation import Conservation
from GUD.ORM.region import Region
from GUD.ORM.source import Source
from . import ParseUtils

usage_msg = """
usage: %s --genome STR --version STR [-h] [options]
""" % os.path.basename(__file__)

help_msg = """%s
inserts features from the UCSC's "phastConsElements100way" table
into GUD.

  --genome STR        genome assembly
  --version STR       e.g. "2015-05-08" (date last updated)

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
# Functions   #
#-------------#

def parse_args():
    """
    This function parses arguments provided via the command line and returns an {argparse} object.
    """

    parser = argparse.ArgumentParser(add_help=False)

    # Mandatory args
    parser.add_argument("--genome")
    parser.add_argument("--version")

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
    if not args.genome or not args.version:
        error = ["%s\n%s" % (usage_msg, os.path.basename(__file__)), "error", "argument \"--genome\" \"--version\" are required\n"]
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
    GUDUtils.pwd  = args.pwd
    GUDUtils.host = args.host
    GUDUtils.port = args.port
    GUDUtils.db = args.db

    # Insert conservation data
    conservation_to_gud(args.genome, args.version, args.dummy_dir, args.remove, args.test, args.threads)

def conservation_to_gud(genome, version, dummy_dir="/tmp/", remove=False, test=False, threads=1):
    """
    e.g. python -m GUD.parsers.phastconselem2gud --genome hg38 --version abcd --test -P 3306
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

    # Create dummy dir
    subdir = "%s.%s" % (genome, os.path.basename(__file__))
    dummy_dir = os.path.join(dummy_dir, subdir)
    if not os.path.isdir(dummy_dir):
        os.makedirs(dummy_dir)

    # Download data
    data_file, url = _download_data(genome, dummy_dir)

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
    ParseUtils.create_table(Conservation)

    # Start a new session
    session = Session()

    # Get valid chromosomes
    chroms = ParseUtils.get_chroms(session)

    # Get source
    source = Source()
    m = re.search("^%s/*(.+).txt.gz$" % dummy_dir, data_file) 
    source.name = m.group(1)
    source.source_metadata = "%s," % version
    source.metadata_descriptor = "version,"
    source.url = url
    ParseUtils.upsert_source(session, source)
    source = ParseUtils.get_source(session, source.name, source.source_metadata, source.metadata_descriptor, url)

    # This is ABSOLUTELY necessary to prevent MySQL from crashing!
    session.close()
    engine.dispose()

    # Split data
    data_files = _split_data(data_file, threads)

    # Parallelize inserts to the database
    ParseUtils.insert_data_files_in_parallel(data_files, partial(_insert_data, test=test), threads)

    # Remove files
    if remove:
        shutil.rmtree(dummy_dir)

    # Remove session
    Session.remove()

def _download_data(genome, dummy_dir="/tmp/"):

    # Python 3+
    if sys.version_info > (3, 0):
        from urllib.request import urlopen, urlretrieve
    # Python 2.7
    else:
        from urllib import urlopen, urlretrieve

    # Initialize
    ftp_files = []
    url = "ftp://hgdownload.soe.ucsc.edu/goldenPath/%s/database" % genome
    regexp = re.compile("(phastConsElements[0-9]+way.txt.gz)")

    # List files in URL
    u = urlopen(url)
    content = u.read().decode("utf-8")
    for file_name in filter(regexp.search, content.split("\n")):
        m = re.search("(phastConsElements[0-9]+way.txt.gz)", file_name)
        n = re.search("phastConsElements([0-9]+)way.txt.gz", m.group(1))
        ftp_files.append((int(n.group(1)), m.group(1)))

    # Sort files
    ftp_files.sort(key=lambda x: x[0], reverse=True)

    # Download data
    data_file = os.path.join(dummy_dir, ftp_files[0][1])
    if not os.path.exists(data_file):
        f = urlretrieve(os.path.join(url, ftp_files[0][1]), data_file)

    return(data_file, os.path.join(url, ftp_files[0][1]))

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
        region.chrom = str(line[1])
        if region.chrom.startswith("chr"):
            region.chrom = region.chrom[3:]
        if region.chrom not in chroms:
            continue
        region.start = int(line[2])
        region.end = int(line[3])
        region.bin = assign_bin(region.start, region.end)
        ParseUtils.upsert_region(session, region)

        # Get region
        region = ParseUtils.get_region(session, region.chrom, region.start, region.end)

        # Upsert conservation
        conservation = Conservation()
        conservation.region_id = region.uid
        conservation.source_id = source.uid
        conservation.score = float(line[5])
        ParseUtils.upsert_conservation(session, conservation)

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