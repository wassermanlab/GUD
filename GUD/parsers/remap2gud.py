#!/usr/bin/env python

import argparse
from binning import assign_bin
from copy import copy
from functools import partial
import getpass
from multiprocessing import cpu_count
import os
from pybedtools import BedTool, cleanup, set_tempdir
import re
import subprocess
import sys
import tarfile

# Import from GUD module
from GUD import GUDUtils
from GUD.ORM.experiment import Experiment
from GUD.ORM.region import Region
from GUD.ORM.sample import Sample
from GUD.ORM.source import Source
from GUD.ORM.tf_binding import TFBinding
from . import ParseUtils

usage_msg = """
usage: %s --genome STR --samples FILE [-h] [options]
""" % os.path.basename(__file__)

help_msg = """%s
inserts ReMap transcription factor ChIP-seq data into GUD.

  --genome STR        genome assembly
  --samples FILE      ReMap samples (manually-curated)

optional arguments:
  -h, --help          show this help message and exit
  --dummy-dir DIR     dummy directory (default = "/tmp/")
  -m, --merge         merge genomic regions using bedtools
                      (default = False)
  -t, --test          limit the total of inserts to ~1K per
                      thread for testing (default = False)
  --threads INT       number of threads to use (default = %s)

mysql arguments:
  -d STR, --db STR    database name (default = "%s")
  -H STR, --host STR  host name (default = "localhost")
  -p STR, --pwd STR   password (default = ignore this option)
  -P STR, --port STR  port number (default = %s)
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
    parser.add_argument("--samples")

    # Optional args
    optional_group = parser.add_argument_group("optional arguments")
    optional_group.add_argument("-h", "--help", action="store_true")
    optional_group.add_argument("--dummy-dir", default="/tmp/")
    optional_group.add_argument("-m", "--merge", action="store_true")
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
    if not args.genome or not args.samples:
        error = ["%s\n%s" % (usage_msg, os.path.basename(__file__)), "error", "arguments \"--genome\" \"--samples\" are required\n"]
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

def main():

    # Parse arguments
    args = parse_args()

    # Set MySQL options
    GUDUtils.user = args.user
    GUDUtils.pwd = args.pwd
    GUDUtils.host = args.host
    GUDUtils.port = args.port
    GUDUtils.db = args.db

    # Insert ReMap data
    remap_to_gud(args.genome, args.samples, args.dummy_dir, args.merge, args.test, args.threads)

def remap_to_gud(genome, samples_file, dummy_dir="/tmp/", merge=False, test=False, threads=1):
    """
    e.g. python -m GUD.parsers.remap2gud --genome hg38 --samples ./samples/ReMap.tsv --dummy-dir ./tmp/ --merge --test -P 3306
    """

    # Initialize
    global chroms
    global engine
    global experiment
    global Session
    global samples
    global source
    experiment_type = "ChIP-seq"
    source_name = "ReMap"
    set_tempdir(dummy_dir) # i.e. for pyBedTools

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
    ParseUtils.create_table(TFBinding)

    # Start a new session
    session = Session()

    # Get valid chromosomes
    chroms = ParseUtils.get_chroms(session)

    # Get experiment
    experiment = Experiment()
    experiment.name = experiment_type
    ParseUtils.upsert_experiment(session, experiment)
    experiment = ParseUtils.get_experiment(session, experiment_type)

    # Get samples
    samples = _get_samples(session, samples_file)

    # Get source
    source = Source()
    source.name = source_name
    ParseUtils.upsert_source(session, source)
    sources = ParseUtils.get_source(session, source_name)
    source = next(iter(sources))

    # This is ABSOLUTELY necessary to prevent MySQL from crashing!
    session.close()
    engine.dispose()

    # Download data
    download_file = _download_data(genome, dummy_dir)

    # For each BED file...
    for bed_file in _unwind_bed_files(download_file, dummy_dir):

        # Prepare data
        data_file = _preprocess_data(bed_file, dummy_dir, merge)

        # Split data
        data_files = _split_data(data_file, threads)

        # Parallelize inserts to the database
        ParseUtils.insert_data_files_in_parallel(data_files, partial(_insert_data, test=test), threads)

        # Remove data files
        if not test:
            if os.path.exists(data_file):
                os.remove(data_file)
            for data_file in data_files:
                if os.path.exists(data_file):
                    os.remove(data_file)

    # Remove downloaded file
    if os.path.exists(download_file) and not test:
        os.remove(download_file)

    # Dispose session
    Session.remove()

def _get_samples(session, file_name):

    # Initialize
    samples = {}
    bugged_samples = {
        "LNCAPABL": "LNCAP"
    }

    # For each line...
    for line in ParseUtils.parse_tsv_file(file_name):

        # If add...
        if line[7] == "Yes":

            # Get sample
            sample = Sample()
            sample.name = line[6]
            sample.treatment = False
            if line[8] == "Yes":
                sample.treatment = True
            sample.cell_line = False
            if line[9] == "Yes":
                sample.cell_line = True
            sample.cancer = False
            if line[10] == "Yes":
                sample.cancer = True

            # Upsert sample
            ParseUtils.upsert_sample(session, sample)

            # Get sample ID
            sample = ParseUtils.get_sample(session, sample.name, sample.X, sample.Y, sample.treatment, sample.cell_line, sample.cancer)

            # Add sample
            samples.setdefault(line[0], sample.uid)

    # For each bugged sample...
    for sample in bugged_samples:
        samples[sample] = samples[bugged_samples[sample]]

    return(samples)

def _download_data(genome, dummy_dir="/tmp/"):

    # Python 3+
    if sys.version_info > (3, 0):
        from urllib.request import urlretrieve
    # Python 2.7
    else:
        from urllib import urlretrieve

    # Initialize
    url = "http://tagc.univ-mrs.fr/remap/download/remap2018/%s/MACS" % genome
    ftp_file = "remap2018_TF_archive_all_macs2_%s_v1_2.tar.gz" % genome

    # Download data
    dummy_file = os.path.join(dummy_dir, ftp_file)
    if not os.path.exists(dummy_file):
        f = urlretrieve(os.path.join(url, ftp_file), dummy_file)

    return(dummy_file)

def _unwind_bed_files(file_name, dummy_dir="/tmp/"):
    """
    Unwinds tarball into 1x BED file per TF.
    """

    with tarfile.open(file_name) as tar:

        # For each member...
        for member in tar.getmembers():

            # Skip if BED file already exists
            bed_file = os.path.join(dummy_dir, member.name)
            if not os.path.exists(bed_file):
                tar.extract(member, dummy_dir)

            yield(bed_file)

def _preprocess_data(file_name, dummy_dir="/tmp/", merge=False):

    # Initialize
    dummy_files = []
    m = re.search("remap2018_(.+)_all_macs2", file_name)
    tf = m.group(1)

    # Skip if intersection file already exists
    bed_file = os.path.join(dummy_dir, "%s.bed" % tf)
    if not os.path.exists(bed_file):

        # Sort BED
        dummy_file = os.path.join(dummy_dir, "dummy.sorted.bed")
        if not os.path.exists(dummy_file):

            # UNIX sort is more efficient than bedtools
            cmd = "zless %s | cut -f 1-4 | sort -T %s -k1,1 -k2,2n > %s" % (file_name, dummy_dir, dummy_file)
            subprocess.call(cmd, shell=True)

        # Add dummy file
        dummy_files.append(dummy_file)

        # Merge BED
        if merge:

            # Initialize
            a = BedTool(dummy_file)

            # Skip if already merged
            dummy_file = os.path.join(dummy_dir, "dummy.merged.bed")
            if not os.path.exists(dummy_file):
                a.merge(stream=True).saveas(dummy_file)

            # Add dummy file
            dummy_files.append(dummy_file)

        # Intersect
        a = BedTool(dummy_files[0])
        b = BedTool(dummy_files[-1])
        a.intersect(b, wa=True, wb=True, stream=True).saveas(bed_file)

        # Clean PyBedTools files
        cleanup(remove_all=True)

        # Remove ALL dummy files
        for dummy_file in dummy_files:
            os.remove(dummy_file)

    return(bed_file)

def _split_data(data_file, threads=1):

    # Initialize
    split_files = []

    # For each chromosome...
    for chrom in chroms:

        # Skip if file already split
        split_file = "%s.%s" % (data_file, chrom)
        if not os.path.exists(split_file):

            # Parallel split
            cmd = 'zless %s | parallel -j %s --pipe --block 2M -k grep "^%s[[:space:]]" > %s' % (data_file, threads, chrom, split_file)
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

        # Skip empty lines
        if not line:
            continue

        # Get region
        region = Region()
        region.chrom = line[-3]
        region.start = int(line[-2])
        region.end = int(line[-1])
        region.bin = assign_bin(region.start, region.end)

        # Ignore non-standard chroms, scaffolds, etc.
        if region.chrom not in chroms:
            continue

        # Get TF, sample
        m = re.search("^\S+\.(\S+)\.(\S+)$", line[3])
        n = re.search("^([^_]+)", m.group(1))
        tf_name = n.group(1)
        n = re.search("^([^_]+)", m.group(2))
        sample_name = n.group(1)

        # Ignore samples
        if sample_name not in samples:
            continue

        # Upsert region
        ParseUtils.upsert_region(session, region)

        # Get region ID
        region = ParseUtils.get_region(session, region.chrom, region.start, region.end, region.strand)

        # Get TF
        tf = TFBinding()
        tf.region_id = region.uid
        tf.sample_id = samples[sample_name]
        tf.experiment_id = experiment.uid
        tf.source_id = source.uid
        tf.tf = tf_name

        # Upsert tf
        ParseUtils.upsert_tf(session, tf)

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
