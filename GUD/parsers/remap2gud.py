#!/usr/bin/env python

import argparse
from binning import assign_bin
from functools import partial
import getpass
from multiprocessing import cpu_count
import os
import pickle
import re
import shutil
import subprocess
import sys
import tarfile

# Import from GUD module
from GUD import GUDUtils
from GUD.ORM.experiment import Experiment
from GUD.ORM.gene import Gene
from GUD.ORM.region import Region
from GUD.ORM.sample import Sample
from GUD.ORM.source import Source
from GUD.ORM.tf_binding import TFBinding
from . import ParseUtils

usage_msg = """
usage: %s --genome STR --samples FILE [-h] [options]
""" % os.path.basename(__file__)

help_msg = """%s
inserts ReMap 2020 transcription factor ChIP-seq data into GUD.

  --genome STR        genome assembly
  --samples FILE      ReMap samples (manually-curated)

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
    parser.add_argument("--samples")

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

    # Insert ReMap data
    remap_to_gud(args.genome, args.samples, args.dummy_dir, args.remove, args.test, args.threads)

def remap_to_gud(genome, samples_file, dummy_dir="/tmp/", remove=False, test=False, threads=1):
    """
    e.g. python -m GUD.parsers.remap2gud --genome hg38 --samples ./samples/ReMap.tsv
    """

    # Initialize
    global chroms
    global engine
    global experiment
    global genes
    global Session
    global samples
    global source
    experiment_type = "ChIP-seq"
    source_name = "ReMap"

    # Testing
    if test:
        global current_process
        from multiprocessing import current_process

    # Create dummy dir
    subdir = "%s.%s" % (genome, os.path.basename(__file__))
    dummy_dir = os.path.join(dummy_dir, subdir)
    if not os.path.isdir(dummy_dir):
        os.makedirs(dummy_dir)

    # Unless pickle file exists
    pickle_file = os.path.join(dummy_dir, "datasets.pickle")
    if not os.path.exists(pickle_file):
        # Get ReMap datasets
        datasets = _get_ReMap_datasets()
        handle = ParseUtils._get_file_handle(pickle_file, "wb")
        pickle.dump(datasets, handle)

    else:
        handle = ParseUtils._get_file_handle(pickle_file, "rb")
        datasets = pickle.load(handle)

    # Group ReMap datasets by target name
    grouped_datasets = _group_ReMap_datasets(datasets)

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
    ParseUtils.create_table(TFBinding)

    # Start a new session
    session = Session()

    # Get valid chromosomes
    chroms = ParseUtils.get_chroms(session)

    # Get all genes
    q = Gene().get_all_gene_symbols(session)
    genes = set([g[0] for g in q.all()])

    # Get experiment
    experiment = Experiment()
    experiment.name = experiment_type
    ParseUtils.upsert_experiment(session, experiment)
    experiment = ParseUtils.get_experiment(session, experiment_type)

    # This is ABSOLUTELY necessary to prevent MySQL from crashing!
    session.close()
    engine.dispose()

    # Get samples
    samples = _get_samples(samples_file)

    # Insert samples/sources
    _insert_samples_and_sources(samples, datasets)

    # # Get source
    # source = Source()
    # source.name = source_name
    # ParseUtils.upsert_source(session, source)
    # sources = ParseUtils.get_source(session, source_name)
    # source = next(iter(sources))



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

def _get_ReMap_datasets():

    # Initialize
    datasets = {}

    import coreapi
    import json

    client = coreapi.Client()
    codec = coreapi.codecs.CoreJSONCodec()

    url = "http://remap.univ-amu.fr:80/api/v1/datasets/findByTaxid/taxid=9606"

    for dataset in json.loads(codec.encode(client.get(url)))["datasets"]:
        datasets.setdefault(dataset["dataset_name"], dataset)

    return(datasets)

def _group_ReMap_datasets(datasets):

    # Initialize
    grouped_datasets = {}

    # For each accession...
    for dataset in datasets:

        # Initialize
        target_name = datasets[dataset]["target_name"]

        # Group datasets
        grouped_datasets.setdefault(target_name, set())
        grouped_datasets[target_name].add(dataset)

    return(grouped_datasets)

# def _download_data(genome, dummy_dir="/tmp/"):

#     # Initialize
#     dummy_files = []

#     # Python 3+
#     if sys.version_info > (3, 0):
#         from urllib.request import urlretrieve
#     # Python 2.7
#     else:
#         from urllib import urlretrieve

#     if genome == "hg19":

#         url = "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/"
#         chains_file = "hg38ToHg19.over.chain.gz"
#         dummy_files.append(os.path.join(dummy_dir, chains_file))

#         # Download data
#         if not os.path.exists(dummy_files[-1]):
#             f = urlretrieve(os.path.join(url, chains_file), dummy_files[-1])

#     else:
#         dummy_files.append(None)

#     # Download data
#     url = "http://remap.univ-amu.fr/storage/remap2020/hg38/MACS2"
#     ftp_file = "remap2020_all_macs2_hg38_v1_0.bed.gz"
#     dummy_files.append(os.path.join(dummy_dir, ftp_file))
#     if not os.path.exists(dummy_files[-1]):
#         f = urlretrieve(os.path.join(url, ftp_file), dummy_files[-1])

#     return(dummy_files[::-1], url)

def _get_samples(file_name):

    # Initialize
    samples = {}

    # For each line...
    for line in ParseUtils.parse_tsv_file(file_name):

        # Initialize
        sample_name = line[0]
        treatment = False
        cell_line = False
        cancer = False
        sex = None
        add = False

        if line[2] == "Yes":
            treatment = True
        if line[3] == "Yes":
            cell_line = True
        if line[4] == "Yes":
            cancer = True
        if line[5] == "female" or line[5] == "male":
            sex = line[5]
        if line[6] == "Yes":
            add = True

        # Add sample
        samples.setdefault(sample_name, (treatment, cell_line, cancer, sex, add))

    return(samples)

def _insert_samples_and_sources(samples, datasets):

    # Initialize
    session = Session()
    i_have_been_warned = False

    # For each dataset...
    for dataset in sorted(datasets):

        # Initialize
        biotype_name = datasets[dataset]["biotype_name"]

        if biotype_name not in samples:
            warnings.warn("missing sample: %s" % biotype_name, Warning, stacklevel=2)
            i_have_been_warned = True
            continue

        if not samples[biotype_name][-1]:
            continue

        # Upsert sample
        sample = Sample()
        sample.name = biotype_name
        sample.treatment = samples[biotype_name][0]
        sample.cell_line = samples[biotype_name][1]
        sample.cancer = samples[biotype_name][2]
        if samples[biotype_name][3] == "female":
            sample.X = 2
            sample.Y = 0
        if samples[biotype_name][3] == "male":
            sample.X = 1
            sample.Y = 1
        ParseUtils.upsert_sample(session, sample)

    if i_have_been_warned:
        error = ["%s" % os.path.basename(__file__), "error", "missing samples"]
        print(": ".join(error))
        exit(0)

    exit(0)
        # # Upsert sample
        # sample = Sample()
        # if not encodes[accession].summary:
        #     sample.name = encodes[accession].biosample_name
        # else:
        #     sample.name = encodes[accession].summary
        # sample.treatment = encodes[accession].treatment
        # sample.cell_line = encodes[accession].cell_line
        # sample.cancer = encodes[accession].cancer
        # if encodes[accession].sex is not None:
        #     sample.X = encodes[accession].X
        #     sample.Y = encodes[accession].Y
        # ParseUtils.upsert_sample(session, sample)

        # # Upsert source
        # source = Source()
        # source.name = "ENCODE"
        # # This is to satisfy ENCODE's citing guidelines:
        # # https://www.encodeproject.org/help/citing-encode/
        # source.source_metadata = "%s," % accession
        # source.metadata_descriptor = "accession,"
        # source.url = encodes[accession].download_url
        # ParseUtils.upsert_source(session, source)

    # This is ABSOLUTELY necessary to prevent MySQL from crashing!
    session.close()

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
        region = ParseUtils.get_region(session, region.chrom, region.start, region.end)

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
