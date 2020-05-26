#!/usr/bin/env python

import argparse
from binning import assign_bin
from functools import partial
import getpass
from multiprocessing import cpu_count
import os
from pybedtools import BedTool
import re
import shutil
import subprocess
import sys
import warnings
from zipfile import ZipFile

# Import from GUD module
from GUD import GUDUtils
from GUD.ORM.dna_accessibility import DNAAccessibility
from GUD.ORM.experiment import Experiment
from GUD.ORM.region import Region
from GUD.ORM.sample import Sample
from GUD.ORM.source import Source
from . import ParseUtils

usage_msg = """
usage: %s --genome STR [-h] [options]
""" % os.path.basename(__file__)

help_msg = """%s
inserts features from the Brain Open Chromatin Atlas (BOCA)
into GUD.

  --genome STR        genome assembly

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
    if not args.genome:
        error = ["%s\n%s" % (usage_msg, os.path.basename(__file__)), "error", "argument \"--genome\" is required\n"]
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

    # Insert BOCA data
    boca_to_gud(args.genome, args.dummy_dir, args.remove, args.test, args.threads)

def boca_to_gud(genome, dummy_dir="/tmp/", remove=False, test=False, threads=1):
    """
    e.g. python -m GUD.parsers.boca2gud --genome hg38 --test -P 3306
    """

    # Initialize
    global chroms
    global engine
    global experiment
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
    ParseUtils.create_table(DNAAccessibility)

    # Start a new session
    session = Session()

    # Get valid chromosomes
    chroms = ParseUtils.get_chroms(session)

    # Get experiment
    experiment = Experiment()
    experiment.name = "ATAC-seq"
    ParseUtils.upsert_experiment(session, experiment)
    experiment = ParseUtils.get_experiment(session, experiment.name)    

    # Get source
    source = Source()
    source.name = "BOCA"
    source.url = url+data_files[0]
    ParseUtils.upsert_source(session, source)
    source = ParseUtils.get_source(session, source.name, source.source_metadata,
                                   source.metadata_descriptor, source.url)

    # Prepare data
    data_file = _preprocess_data(session, data_files, dummy_dir)

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

    # Initialize
    data_files = []
    ftp_files = []

    # Python 3+
    if sys.version_info > (3, 0):
        from urllib.request import urlretrieve
    # Python 2.7
    else:
        from urllib import urlretrieve

    if genome == "hg19":
        ftp_files.append("boca_peaks.zip")
        ftp_files.append("boca_peaks_consensus_no_blacklisted_regions.bed")
    else:
        ftp_files.append("boca_peaks_hg38_lifted.zip")
        ftp_files.append("boca_peaks_consensus_no_blacklisted_regions_hg38_lifted.bed")

    # Download data
    url = "https://bendlj01.u.hpc.mssm.edu/multireg/resources/"
    for data_file in ftp_files:
        data_files.append(os.path.join(dummy_dir, data_file)) 
        if not os.path.exists(data_files[-1]):
            f = urlretrieve(os.path.join(url, data_file), data_files[-1])

    return(data_files, url)

def _preprocess_data(session, data_files, dummy_dir="/tmp/"):

    # Initialize
    abbreviations = _get_abbreviations()


    # Skip if BED file exists
    bed_file = os.path.join(dummy_dir, "BOCA.bed")
    if not os.path.exists(bed_file):

        # Initialize
        lines = []

        # For each zip file...
        with ZipFile(data_files[0], "r") as z:
            for zip_info in z.filelist:
                m = re.search("(\w+)\_(glia|neuron).bed$", zip_info.filename)
                if not m:
                    continue
                if m.group(1) not in abbreviations:
                    continue
                s = "%s_%s" % (m.group(1), m.group(2))
                for line in z.read(zip_info.filename).splitlines():
                    chrom, start, end = line.decode().split("\t")[:3]
                    lines.append("%s\t%s\t%s\t%s" % (chrom, start, end, s))
                # Upsert sample
                sample = Sample()
                sample.name = "%s (%s)" % (abbreviations[m.group(1)], m.group(2))
                sample.treatment = False
                sample.cell_line = False
                sample.cancer = False
                ParseUtils.upsert_sample(session, sample)

        # Get intersections
        a = BedTool("\n".join(lines), from_string=True).sort()
        b = BedTool(data_files[1]).sort()
        a = a.intersect(b, wa=True, u=True, sorted=True)

        # Save
        a.saveas(bed_file)

    return(bed_file)

def _get_abbreviations():

    # Initialize
    abbreviations = {
        "DLPFC": "dorsolateral prefrontal cortex",
        "OFC": "orbitofrontal cortex",
        "VLPFC": "ventrolateral prefrontal cortex",
        "ACC": "anterior cingulate cortex",
        "STC": "superior temporal cortex",
        "ITC": "inferior temporal cortex",
    	"PMC": "primary motor cortex",
        "INS": "insula",
        "PVC": "primary visual cortex",
        "AMY": "amygdala",
        "HIP": "hippocampus",
        "MDT": "mediodorsal thalamus",
        "NAC": "nucleus accumbens",
        "PUT": "putamen"
    }

    return(abbreviations)

def _split_data(data_file, threads=1):

    # Initialize
    split_files = []
    split_dir = os.path.dirname(os.path.realpath(data_file))

    # Get number of lines
    output = subprocess.check_output(["wc -l %s" % data_file], shell=True)
    m = re.search("(\d+)", str(output))
    L = float(m.group(1))

    # Split
    prefix = "%s." % data_file.split("/")[-1]
    cmd = "split -d -l %s %s %s" % (int(L/threads)+1, data_file, os.path.join(split_dir, prefix))
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
    samples = {}
    abbreviations = _get_abbreviations()

    # Testing
    if test:
        lines = 0
        print(current_process().name)

    # For each line...
    for line in ParseUtils.parse_tsv_file(data_file):

        # Initialize
        abbrev, origin = line[3].split("_")

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

        # Get sample
        s = "%s (%s)" % (abbreviations[abbrev], origin)
        if s not in samples:
            sample = Sample()
            sample.name = s
            sample.treatment = False
            sample.cell_line = False
            sample.cancer = False
            sample = ParseUtils.get_sample(session, sample.name, sample.X, sample.Y, sample.treatment, sample.cell_line, sample.cancer)
            samples.setdefault(s, sample)
        else:
            sample = samples[s]

        # Upsert DNA accessibility
        dna_accessibility = DNAAccessibility()
        dna_accessibility.region_id = region.uid
        dna_accessibility.sample_id = sample.uid
        dna_accessibility.experiment_id = experiment.uid
        dna_accessibility.source_id = source.uid
        dna_accessibility.score = None
        dna_accessibility.peak = None
        ParseUtils.upsert_accessibility(session, dna_accessibility)

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
