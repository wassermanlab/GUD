#!/usr/bin/env python

import argparse
from binning import assign_bin
import getpass
from multiprocessing import cpu_count
from multiprocessing import current_process
import os
from pybedtools import BedTool, cleanup, set_tempdir
import re
from sqlalchemy_utils import database_exists
import sys
import tarfile
# Python 3+
if sys.version_info > (3, 0):
    from urllib.request import urlretrieve
# Python 2.7
else:
    from urllib import urlretrieve

# Import from GUD module
from GUD import GUDglobals
from GUD.ORM.experiment import Experiment
from GUD.ORM.region import Region
from GUD.ORM.sample import Sample
from GUD.ORM.source import Source
from GUD.ORM.tf_binding import TFBinding
from . import _get_chroms, _get_db_name, _get_experiment, _get_region, _get_sample, _get_source, _initialize_gud_db, _initialize_engine_session, _process_data_in_chunks, _upsert_experiment, _upsert_region, _upsert_sample, _upsert_source, _upsert_tf

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
    parser.add_argument("--samples")

    # Optional args
    optional_group = parser.add_argument_group("optional arguments")
    optional_group.add_argument("-h", "--help", action="store_true")
    optional_group.add_argument("--dummy-dir", default="/tmp/")
    optional_group.add_argument("-m", "--merge", action="store_true")
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
    if not args.genome or not args.samples:
        error = ["%s\n%s" % (usage_msg, os.path.basename(__file__)), "error", "arguments \"--genome\" \"--samples\" are required\n"]
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

    # Insert ReMap data
    remap_to_gud(args.user, args.pwd, args.host, args.port, args.db, args.genome, args.samples, args.merge, args.dummy_dir, args.threads)

def remap_to_gud(user, pwd, host, port, db, genome, samples_file, merge=False, dummy_dir="/tmp/", threads=1):
    """
    python -m GUD.parsers.remap2gud --genome hg19 --samples ./samples/ReMap.tsv --dummy-dir ./tmp/ --merge
    """

    # Globals
    global chroms
    global engine
    global experiment
    global Session
    global samples
    global source

    # Initialize
    experiment_type = "ChIP-seq"
    source_name = "ReMap"
    set_tempdir(dummy_dir) # i.e. for pyBedTools

    # Download data
    dummy_file = _download_data(genome, dummy_dir)

    # Get database name
    db_name = _get_db_name(user, pwd, host, port, db)

    # If database does not exist...
    if not database_exists(db_name):
        _initialize_gud_db(user, pwd, host, port, db, genome)

    # Get engine/session
    engine, Session = _initialize_engine_session(db_name)
    session = Session()

    # Create table
    if not engine.has_table(TFBinding.__tablename__):
        TFBinding.__table__.create(bind=engine)

    # Get valid chromosomes
    chroms = _get_chroms(session)

    # Get experiment
    experiment = Experiment()
    experiment.name = experiment_type
    _upsert_experiment(session, experiment)
    experiment = _get_experiment(session, experiment_type)    

    # Get samples
    samples = _get_samples(session, samples_file)

    # Get source
    source = Source()
    source.name = source_name
    _upsert_source(session, source)
    source = _get_source(session, source_name)

    # This is ABSOLUTELY necessary to prevent MySQL from crashing!
    session.close()
    engine.dispose()

    # For each BED file...
    for bed_file in _unwind_bed_files(dummy_file, dummy_dir):

        # Prepare data
        data_file = _preprocess_data(bed_file, dummy_dir, merge)

        # Parallelize inserts to the database
        _process_data_in_chunks(data_file, _insert_data_in_chunks, threads)

        # Remove data file
        if os.path.exists(data_file):
            os.remove(data_file)

    # Dispose session
    Session.remove()

    # # Remove downloaded file
    # if os.path.exists(dummy_file):
    #     os.remove(dummy_file)

def _get_samples(session, file_name):

    # Initialize
    samples = {}
    bugged_samples = {
        "LNCAPABL": "LNCAP"
    }

    # For each line...
    for line in GUDglobals.parse_tsv_file(file_name):

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
            _upsert_sample(session, sample)

            # Get sample ID
            sample = _get_sample(session, sample.name, sample.treatment, sample.cell_line, sample.cancer)

            # Add sample
            samples.setdefault(line[0], sample.uid)

    # For each sample...
    for sample in bugged_samples:

        # Fix sample
        samples[sample] = samples[bugged_samples[sample]]

    return(samples)

def _download_data(genome, dummy_dir="/tmp/"):

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
    m = re.search("remap2018_(.+)_all_macs2", file_name)
    tf = m.group(1)

    # Skip if intersection file already exists
    bed_file = os.path.join(dummy_dir, "%s.bed" % tf)
    if not os.path.exists(bed_file):

        # Sort BED
        a = BedTool(file_name)
        a = a.sort()

        # Merge BED
        if merge:
            b = a.merge()
        else:
            b = a

        # Intersect
        a.intersect(b, wa=True, wb=True).saveas(bed_file)

        # Clean PyBedTools files
        cleanup(remove_all=True)

    return(bed_file)

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
        region.chrom = line[9]
        region.start = int(line[10])
        region.end = int(line[11])
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
        _upsert_region(session, region)

        # Get region ID
        region = _get_region(session, region.chrom, region.start, region.end, region.strand)

        # Get TF
        tf = TFBinding()
        tf.regionID = region.uid
        tf.sampleID = samples[sample_name]
        tf.experimentID = experiment.uid
        tf.sourceID = source.uid
        tf.tf = tf_name

        # Upsert tf
        _upsert_tf(session, tf)

    # This is ABSOLUTELY necessary to prevent MySQL from crashing!
    session.close()
    engine.dispose()

#-------------#
# Main        #
#-------------#

if __name__ == "__main__": main()
