#!/usr/bin/env python

import argparse
from binning import assign_bin
from functools import partial
import getpass
from multiprocessing import Pool, cpu_count
import os
import pickle
import re
import shutil
import subprocess
import sys
import tarfile
# Python 3+
if sys.version_info > (3, 0):
    from urllib.request import urlretrieve
# Python 2.7
else:
    from urllib import urlretrieve

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
    global dataset2samplesNsourcesNtfs
    global engine
    global experiment
    global genes
    global Session
    global samples
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

    # LiftOver
    if genome == "hg19":
        liftOver = True
        chains_file = _get_chains_file(dummy_dir)
    else:
        liftOver = False
        chains_file = None

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
    experiment.name = "ChIP-seq"
    ParseUtils.upsert_experiment(session, experiment)
    experiment = ParseUtils.get_experiment(session, "ChIP-seq")

    # This is ABSOLUTELY necessary to prevent MySQL from crashing!
    session.close()
    engine.dispose()

    # Get samples
    samples = _get_samples(samples_file)

    # Insert samples/sources
    dataset2samplesNsourcesNtfs = _insert_samples_and_sources_and_tfs(samples, datasets, liftOver)
    exit(0)

    # Group ReMap datasets by target name
    grouped_datasets = _group_ReMap_datasets(datasets)

    # Partial function to enable parallelization
    partial_func = partial(_insert_data, test=test, chains_file=chains_file)

    # For each TF...
    for tf in sorted(grouped_datasets):

        # Prepare data
        data_file = _preprocess_data(tf, grouped_datasets, datasets, dummy_dir,
            test, threads)

        # Skip
        if not os.path.exists(data_file):
            continue

        # Split data
        data_files = _split_data(data_file, threads)

        # Parallelize inserts to the database
        ParseUtils.insert_data_files_in_parallel(data_files, partial_func, threads)

    # Remove files
    if remove:
        shutil.rmtree(dummy_dir)

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

def _get_chains_file(dummy_dir="/tmp/"):

    # Initialize
    url = "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/"
    chains_file = "hg38ToHg19.over.chain.gz"
    dummy_file = os.path.join(dummy_dir, chains_file)

    # Download data
    if not os.path.exists(dummy_file):
        f = urlretrieve(os.path.join(url, chains_file), dummy_file)

    return(dummy_file)

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

def _insert_samples_and_sources_and_tfs(samples, datasets, liftOver=False):

    # Initialize
    dataset2samplesNsourcesNtfs = {}
    i_have_been_warned = False
    session = Session()
    source_name = "ReMap"

    # For each dataset...
    for dataset in sorted(datasets):

        # Initialize
        biotype_name = datasets[dataset]["biotype_name"]
        dataset_name = datasets[dataset]["dataset_name"]
        target_name = datasets[dataset]["target_name"]

        # Skip if non-valid sample
        if biotype_name not in samples:
            warnings.warn("missing sample: %s" % biotype_name, Warning, stacklevel=2)
            i_have_been_warned = True
            continue
        if not samples[biotype_name][-1]:
            print(biotype_name)
            continue

        # Skip if non-valid target
        if target_name not in genes:
            continue

        # Skip if target has been modified
        if datasets[dataset]["target_modification"] != "":
            continue

        # Upsert sample
        sample = Sample()
        sample.name = biotype_name
        sample.treatment = samples[biotype_name][0]
        if datasets[dataset]["biotype_modification"] != "":
            sample.treatment = True
        sample.cell_line = samples[biotype_name][1]
        sample.cancer = samples[biotype_name][2]
        if samples[biotype_name][3] == "female":
            sample.X = 2
            sample.Y = 0
        if samples[biotype_name][3] == "male":
            sample.X = 1
            sample.Y = 1
        ParseUtils.upsert_sample(session, sample)
        sample = ParseUtils.get_sample(session, sample.name, sample.X, sample.Y, sample.treatment, sample.cell_line, sample.cancer)
        dataset2samplesNsourcesNtfs.setdefault(dataset_name, [])
        dataset2samplesNsourcesNtfs[dataset_name].append(sample)

        # Upsert source
        source = Source()
        source.name = source_name
        source.source_metadata = "%s,%s," % (dataset_name, liftOver)
        source.metadata_descriptor = "dataset,liftOver,"
        source.url = datasets[dataset]["bed_url"]
        ParseUtils.upsert_source(session, source)
        source = ParseUtils.get_source(session, source.name, source.source_metadata, source.metadata_descriptor, source.url)
        dataset2samplesNsourcesNtfs[dataset_name].append(source)

        # Add TF
        dataset2samplesNsourcesNtfs[dataset_name].append(target_name)

    if i_have_been_warned:
        error = ["%s" % os.path.basename(__file__), "error", "missing samples"]
        print(": ".join(error))
        exit(0)

    # This is ABSOLUTELY necessary to prevent MySQL from crashing!
    session.close()
    engine.dispose()

    return(dataset2samplesNsourcesNtfs)

def _group_ReMap_datasets(datasets):

    # Initialize
    grouped_datasets = {}

    # For each accession...
    for dataset in datasets:

        # Initialize
        target_name = datasets[dataset]["target_name"]

        # Skip if non-valid target
        if target_name not in genes:
            continue

        # Skip if target has been modified
        if datasets[dataset]["target_modification"] != "":
            continue

        # Group datasets
        grouped_datasets.setdefault(target_name, set())
        grouped_datasets[target_name].add(dataset)

    return(grouped_datasets)

def _preprocess_data(tf, grouped_datasets, datasets, dummy_dir="/tmp/",
    test=False, threads=1):

    # Initialize
    dummy_files = []
    label = "ChIP-seq.%s" % tf
    tf_datasets = []

    # Skip if BED file exists
    bed_file = os.path.join(dummy_dir, "%s.bed" % label)
    if not os.path.exists(bed_file):

        # Skip if BED file exists
        dummy_file = os.path.join(dummy_dir, "dummy.bed")
        if not os.path.exists(dummy_file):

            # Get TF datasets
            for dataset in grouped_datasets[tf]:
                tf_datasets.append(datasets[dataset])

            # Get ENCODE BED files
            pool = Pool(processes=threads)
            partial_function = partial(_download_ReMap_bed_file, dummy_dir=dummy_dir, test=test)
            for download_file in pool.imap(partial_function, tf_datasets):

                # Skip
                if download_file is None:
                    continue

                # Initialize
                m = re.search("%s/(\S+).bed.gz$" % dummy_dir, download_file)
                dataset = m.group(1)

                # Concatenate
                cmd = "zless %s | awk '{print $1\"\t\"$2\"\t\"$3\"\t%s\t\"$7\"\t\"$10}' >> %s" % (download_file, dataset, dummy_file)
                subprocess.call(cmd, shell=True)

                # Remove downloaded file
                os.remove(download_file)

            pool.close()

        if os.path.exists(dummy_file):

            # Add dummy file
            dummy_files.append(dummy_file)

            # Sort BED
            dummy_file = os.path.join(dummy_dir, "dummy.sorted.bed")
            if not os.path.exists(dummy_file):

                # UNIX parallel sort is more efficient than bedtools
                cmd = "LC_ALL=C sort --parallel=%s -T %s -k1,1 -k2,2n %s > %s" % (str(threads), dummy_dir, dummy_files[0], dummy_file)
                subprocess.call(cmd, shell=True)

            # Add dummy file
            dummy_files.append(dummy_file)

            # Copy file
            shutil.copy(dummy_files[1], bed_file)

            # Remove ALL dummy files
            for dummy_file in dummy_files:
                os.remove(dummy_file)

    return(bed_file)

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

def _download_ReMap_bed_file(dataset, dummy_dir="/tmp/", test=False):

    try:

        # Initialize
        download_file = os.path.join(dummy_dir, dataset["dataset_name"])

        # Testing
        if test:
            print(current_process().name)

        # Download BED file
        download_file += ".bed.gz"
        if not os.path.exists(download_file):
            urlretrieve(dataset["bed_url"], download_file)

        return(download_file)

    except:
        return(None)

def _insert_data(data_file, test=False, chains_file=None):

    # Initialize
    session = Session()
    if chains_file is not None:
        from pyliftover import LiftOver
        lo = LiftOver(chains_file)

    # Testing
    if test:
        lines = 0
        print(current_process().name)

    # For each line...
    for line in ParseUtils.parse_tsv_file(data_file):

        # Initialize
        dataset_name = line[3]

        # Upsert region
        region = Region()
        region.chrom = str(line[0])
        if region.chrom.startswith("chr"):
            region.chrom = region.chrom[3:]
        if region.chrom not in chroms:
            continue
        region.start = int(line[1])
        region.end = int(line[2])
        if chains_file:
            try:
                chrom = "chr%s" % region.chrom
                start = lo.convert_coordinate(chrom, region.start) # i.e. already 0-based
                region.start = start[0][1]
                end = lo.convert_coordinate(chrom, region.end - 1) # i.e. requires 0-based
                region.end = end[0][1] + 1
            except:
                msg = "position could not be found in new assembly"
                warnings.warn("%s: %s" % (msg, line[1]), Warning, stacklevel=2)
                continue
        region.bin = assign_bin(region.start, region.end)
        ParseUtils.upsert_region(session, region)

        # Get region
        region = ParseUtils.get_region(session, region.chrom, region.start, region.end)

        # Get sample
        sample = dataset2samplesNsourcesNtfs[dataset_name][0]

        # Get source
        source = dataset2samplesNsourcesNtfs[dataset_name][1]

        # Upsert tf
        tf = TFBinding()
        tf.region_id = region.uid
        tf.sample_id = sample.uid
        tf.experiment_id = experiment.uid
        tf.source_id = source.uid
        tf.score = float(line[4])
        tf.peak = int(line[5])
        tf.tf = dataset2samplesNsourcesNtfs[dataset_name][2]
        ParseUtils.upsert_tf(session, tf)

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

# def _unwind_bed_files(file_name, dummy_dir="/tmp/"):
#     """
#     Unwinds tarball into 1x BED file per TF.
#     """

#     with tarfile.open(file_name) as tar:

#         # For each member...
#         for member in tar.getmembers():

#             # Skip if BED file already exists
#             bed_file = os.path.join(dummy_dir, member.name)
#             if not os.path.exists(bed_file):
#                 tar.extract(member, dummy_dir)

#             yield(bed_file)
