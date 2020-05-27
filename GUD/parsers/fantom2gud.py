#!/usr/bin/env python

import argparse
from binning import assign_bin
from functools import partial
import getpass
from multiprocessing import Pool, cpu_count
import os
import re
import shutil
import subprocess
import sys
import warnings

# Import from GUD module
from GUD import GUDUtils
from GUD.ORM.enhancer import Enhancer
from GUD.ORM.experiment import Experiment
from GUD.ORM.region import Region
from GUD.ORM.sample import Sample
from GUD.ORM.source import Source
from GUD.ORM.tss import TSS
from . import ParseUtils

usage_msg = """
usage: %s --genome STR --samples FILE --feature STR [-h] [options]
""" % os.path.basename(__file__)

help_msg = """%s
inserts features from FANTOM5 (Functional Annotation of the
Mammalian Genome 5) into GUD.

  --genome STR        genome assembly
  --samples FILE      FANTOM5 samples (manually-curated)
  --feature STR       type of genomic feature ("enhancer"
                      or "tss")

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
    parser.add_argument("--feature")

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
    if not args.genome or not args.samples or not args.feature:
        error = ["%s\n%s" % (usage_msg, os.path.basename(__file__)), "error", "argument \"--genome\" \"--samples\" \"--feature\" are required\n"]
        print(": ".join(error))
        exit(0)

    # Check for invalid feature
    valid_features = ["enhancer", "tss"]
    if args.feature not in valid_features:
        error = ["%s\n%s" % (usage_msg, os.path.basename(__file__)), "error", "argument \"--feature\"", "invalid choice", "\"%s\" (choose from" % args.feature, "%s)\n" % " ".join(["\"%s\"" % i for i in valid_features])]
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

    # Insert FANTOM5 data
    fantom_to_gud(args.genome, args.samples, args.feature, args.dummy_dir, args.remove, args.test, args.threads)

def fantom_to_gud(genome, samples_file, feat_type, dummy_dir="/tmp/", remove=False, test=False, threads=1):
    """
    e.g. python -m GUD.parsers.fantom2gud --genome hg38 --samples ./samples/FANTOM5.tsv --feature tss
    """

    # Initialize
    global chroms
    global engine
    global experiment
    global samples
    global source
    global Feature
    global Session
    global idx

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
    data_files, url = _download_data(genome, feat_type, dummy_dir)

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
    if feat_type == "enhancer":
        Feature = Enhancer
    else:
        Feature = TSS
    ParseUtils.create_table(Feature)
    if Feature.__tablename__ == "transcription_start_sites":
        global Feature2
        from GUD.ORM.expression import Expression
        ParseUtils.create_table(Expression)
        Feature2 = Expression

    # Start a new session
    session = Session()

    # Get valid chromosomes
    chroms = ParseUtils.get_chroms(session)

    # Get all genes
    if Feature.__tablename__ == "transcription_start_sites":
        global genes
        from GUD.ORM.gene import Gene
        q = Gene().get_all_gene_symbols(session)
        genes = set([g[0] for g in q.all()])

    # Get experiment
    experiment = Experiment()
    experiment.name = "CAGE"
    ParseUtils.upsert_experiment(session, experiment)
    experiment = ParseUtils.get_experiment(session, experiment.name)    

    # Get source
    source = Source()
    source.name = "FANTOM5"
    if genome == "hg19" or genome == "mm9":
        source.source_metadata = "True,"
    else:
        source.source_metadata = "False,"
    source.metadata_descriptor = "liftOver,"
    source.url = url
    ParseUtils.upsert_source(session, source)
    source = ParseUtils.get_source(session, source.name, source.source_metadata,
                                   source.metadata_descriptor, url)

    # Get samples
    samples = _get_samples(session, samples_file)

    # This is ABSOLUTELY necessary to prevent MySQL from crashing!
    session.close()
    engine.dispose()

    # Prepare data
    bed_file, idx = _preprocess_data(data_files, feat_type, dummy_dir, test, threads)

    # Split data
    data_files = _split_data(bed_file, threads)

    # Parallelize inserts to the database
    ParseUtils.insert_data_files_in_parallel(data_files, partial(_insert_data, test=test), threads)

    # Remove files
    if remove:
        shutil.rmtree(dummy_dir)

    # Remove session
    Session.remove()

def _download_data(genome, feat_type, dummy_dir="/tmp/"):

    # Initialize
    dummy_files = []
    ftp_files = []

    # Python 3+
    if sys.version_info > (3, 0):
        from urllib.request import urlretrieve
        from urllib.parse import unquote
    # Python 2.7
    else:
        from urllib import urlretrieve
        from urllib2 import unquote

    if genome == "hg19" or genome == "mm9":

        if genome == "hg19":
            url = "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/"
            chains_file = "hg38ToHg19.over.chain.gz"
            genome = "hg38"
        else:
            url = "http://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/"
            chains_file = "mm10ToMm9.over.chain.gz"
            genome = "mm10"

        dummy_files.append(os.path.join(dummy_dir, chains_file))

        # Download data
        if not os.path.exists(dummy_files[-1]):
            f = urlretrieve(os.path.join(url, chains_file), dummy_files[-1])

    else:
        dummy_files.append(None)

    url = "http://fantom.gsc.riken.jp/5/datafiles/reprocessed/%s_latest/extra/" % genome
    if feat_type == "enhancer":
        url += "enhancer"
        ftp_files.append("F5.%s.enhancers.bed.gz" % genome)
        ftp_files.append("F5.%s.enhancers.expression.usage.matrix.gz" % genome)
    else:
        ftp_files.append("%s_fair+new_CAGE_peaks_phase1and2.bed.gz" % genome) 
        dummy_file = os.path.join(dummy_dir, ftp_files[-1])
        # Download data
        if not os.path.exists(dummy_file):
            urlretrieve(os.path.join(url, "CAGE_peaks", ftp_files[-1]), dummy_file)
        url += "CAGE_peaks_expression"
        ftp_files.append("%s_fair+new_CAGE_peaks_phase1and2_tpm_ann.osc.txt.gz" % genome)

    # Download data
    for ftp_file in ftp_files:
        dummy_files.append(os.path.join(dummy_dir, ftp_file))
        if not os.path.exists(dummy_files[-1]):
            urlretrieve(os.path.join(url, ftp_file), dummy_files[-1])

    return(dummy_files[::-1], os.path.join(url, ftp_files[0]))

def _get_samples(session, file_name):

    # Initialize
    samples = {}

    # For each line...
    for line in ParseUtils.parse_tsv_file(file_name):

        # If add...
        if line[6] == "Yes":

            # Upsert sample
            sample = Sample()
            sample_id = line[0]
            sample.name = line[1]
            if line[2] == "Yes":
                sample.treatment = True
            else:
                sample.treatment = False
            if line[3] == "Yes":
                sample.cell_line = True
            else:
                sample.cell_line = False
            if line[4] == "Yes":
                sample.cancer = True
            else:
                sample.cancer = False
            if line[5] == "female":
                sample.X = 2
                sample.Y = 0
            if line[5] == "male":
                sample.X = 1
                sample.Y = 1
            ParseUtils.upsert_sample(session, sample)

            # Get sample
            sample = ParseUtils.get_sample(session, sample.name, sample.X, sample.Y, sample.treatment, sample.cell_line, sample.cancer)

            # Add sample
            samples.setdefault(sample_id, sample)

    return(samples)

def _preprocess_data(data_files, feat_type, dummy_dir="/tmp/", test=False, threads=1):

    # Initialize
    coords = {}
    dummy_files = []
    if feat_type == "enhancer":
        start_idx = 1
    else:
        start_idx = 7

    # If chains file...
    if data_files[-1] is not None:
        from pyliftover import LiftOver
        lo = LiftOver(data_files[-1])        

    # For each line...
    for line in ParseUtils.parse_tsv_file(data_files[1]):
        chrom, start, end = line[:3]
        if data_files[-1] is not None:
            try:
                l = lo.convert_coordinate(chrom, start) # i.e. already 0-based
                start = l[0][1]
                l = lo.convert_coordinate(chrom, end - 1)  # i.e. requires 0-based
                end = l[0][1] + 1
            except:
                msg = "position could not be found in new assembly"
                warnings.warn("%s: %s:%s-%s" % (msg, chrom, start, end), Warning, stacklevel=2)
                continue
        coords.setdefault(line[3], [chrom, start, end, line[5]])

    # Enhancer file is weird: rows have different number of cols;
    # Read as a normal file and split by tabs
    for line in ParseUtils.parse_file(data_files[0]):

        if feat_type == "enhancer":
            idx = re.findall("CNhs\d{5}", str(line))
            break
        else:
            if line.startswith("00Annotation"):
                idx = re.findall("CNhs\d{5}", str(line))
                break

    # Skip if BED file exists
    bed_file = os.path.join(dummy_dir, "CAGE.%s.bed" % feat_type)
    if not os.path.exists(bed_file):

        # Skip if BED file exists
        dummy_file = os.path.join(dummy_dir, "dummy.bed")
        if not os.path.exists(dummy_file):

            # Enhancer file is weird: rows have different number of cols;
            # Read as a normal file and split by tabs
            for line in ParseUtils.parse_file(data_files[0]):

                line = line.split("\t")

                if line[0] not in coords:
                    continue

                else:
                    if feat_type == "enhancer":
                        txt = "\t".join(map(str, coords[line[0]]+line[start_idx:]))
                    else:
                        txt = "\t".join(map(str, coords[line[0]]+[line[1]]+line[start_idx:]))
                    ParseUtils.write(dummy_file, txt)

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

    return(bed_file, idx)

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
        region.chrom = str(line.pop(0))
        if region.chrom.startswith("chr"):
            region.chrom = region.chrom[3:]
        if region.chrom not in chroms:
            continue
        region.start = int(line.pop(0))
        region.end = int(line.pop(0))
        region.bin = assign_bin(region.start, region.end)
        ParseUtils.upsert_region(session, region)

        # Get region
        region = ParseUtils.get_region(session, region.chrom, region.start, region.end)

        # Upsert/get feature
        if Feature.__tablename__ == "transcription_start_sites":
            sample_ids = []
            expression_levels = []
            feature = Feature()
            feature.region_id = region.uid
            feature.experiment_id = experiment.uid
            feature.source_id = source.uid
            feature.strand = line.pop(0)
            try:
                m = re.search("p(\d+)@(\w+)", line.pop(0))
                if m.group(2) in genes:
                    feature.tss = int(m.group(1))
                    feature.gene = m.group(2)
            except:
                feature.tss = 1
                feature.gene = None
            for i in range(len(idx)):
                if float(line[i]) > 0 and idx[i] in samples:
                    sample_ids.append(samples[idx[i]].uid)
                    expression_levels.append(float("{:.3f}".format(line[i])))
            sample_id = "%s," % ",".join(map(str, sample_ids))
            feature.sample_id = sample_id.encode(encoding="UTF-8")
            expression_level = "%s," % ",".join(map(str, expression_levels))
            feature.expression_level = expression_level.encode(encoding="UTF-8")
            ParseUtils.upsert_tss(session, feature)
            feature = ParseUtils.get_tss(session, feature.region_id, feature.experiment_id,
                                         feature.source_id, feature.gene, feature.tss)
            # Upsert feature2
            if feature.gene:
                for i in range(len(sample_ids)):
                    feature2 = Feature2()
                    feature2.expression_level = expression_levels[i]
                    feature2.tss_id = feature.uid
                    feature2.sample_id = sample_ids[i]
                    ParseUtils.upsert_expression(session, feature2)
        else:
            strand = line.pop(0)
            for i in range(len(idx)):
                if int(line[i]) > 0 and idx[i] in samples:
                    feature = Feature()
                    feature.region_id = region.uid
                    feature.experiment_id = experiment.uid
                    feature.source_id = source.uid
                    feature.sample_id = samples[idx[i]].uid
                    ParseUtils.upsert_enhancer(session, feature)

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
