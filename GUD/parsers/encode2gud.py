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
# Python 3+
if sys.version_info > (3, 0):
    from urllib.request import urlretrieve
# Python 2.7
else:
    from urllib import urlretrieve

# Import from GUD module
from GUD import GUDUtils
from GUD.ORM.dna_accessibility import DNAAccessibility
from GUD.ORM.experiment import Experiment
from GUD.ORM.histone_modification import HistoneModification
from GUD.ORM.region import Region
from GUD.ORM.sample import Sample
from GUD.ORM.source import Source
from GUD.ORM.tf_binding import TFBinding
from . import ParseUtils

usage_msg = """
usage: %s --genome STR [-h] [options]
""" % os.path.basename(__file__)

help_msg = """%s
inserts features from ENCODE (Encyclopedia of DNA Elements)
into GUD.

  --genome STR        genome assembly
  --samples FILE      ENCODE samples (manually-curated)
  --feature STR       type of genomic feature ("atac-seq"
                      "accessibility", "histone" or "tf")

optional arguments:
  -h, --help          show this help message and exit
  --dummy-dir DIR     dummy directory (default = "/tmp/")
  -m, --merge         merge genomic regions using bedtools
                      (default = False)
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
# Classes     #
#-------------#

class EncodeMetadata:

    def __init__(self, accession, biosample, download_url, experiment_accession, experiment_type, experiment_target, genome_assembly, output_format, output_type, status, treatments):
        """
        m = re.search(
            "^(3xFLAG|eGFP)?-?(.+)-(human|mouse)$",
            
        )
        if m:
            tag = m.group(1)
            experiment_target = m.group(2)
        """

        # Fix hg38
        if genome_assembly == "GRCh38":
            genome_assembly = "hg38"

        self.accession = accession
        self.biosample = biosample
        self.download_url = download_url
        self.experiment_accession = experiment_accession
        self.experiment_type = experiment_type
        self.experiment_target = experiment_target
        self.genome_assembly = genome_assembly
        self.output_format = output_format
        self.output_type = output_type
        self.status = status
        self.treatments = treatments

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
    optional_group.add_argument("-m", "--merge", action="store_true")
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
    valid_features = ["atac-seq", "accessibility", "histone", "tf"]
    if args.feature not in valid_features:
        error = ["%s\n%s" % (usage_msg, os.path.basename(__file__)), "error", "argument \"--feature\"", "invalid choice", "\"%s\" (choose from" % args.feature, "%s)\n" % " ".join(["\"%s\"" % i for i in valid_features])]
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

    # Insert ENCODE data
    encode_to_gud(args.genome, args.samples, args.feature, args.merge, args.dummy_dir, args.threads)

def encode_to_gud(genome, samples_file, feat_type, merge=False, dummy_dir="/tmp/", threads=1):
    """
    python -m GUD.parsers.encode2gud --genome hg19 --samples --dummy-dir ./tmp/
    """

    # Globals
    global chroms
    global engine
    global experiment
    global Session
    global samples
    global source
    global table

    # Initialize
    source_name = "ENCODE"
    set_tempdir(dummy_dir) # i.e. for pyBedTools

    # Download metadata
    metadata_file = _download_metadata(genome, feat_type, dummy_dir)

    # Get database name
    db_name = GUDUtils._get_db_name()

    # If database does not exist...
    if not database_exists(db_name):
        ParseUtils.initialize_gud_db(db_name, genome)

    # Get engine/session
    engine, Session = GUDUtils._get_engine_session(db_name)
    session = Session()

    # Create table
    if feat_type == "accessibility":
        table = DNAAccessibility()
    elif feat_type == "histone":
        table = HistoneModification()
    else:
        table = TFBinding()
    if not engine.has_table(table.__tablename__):
        table.__table__.create(bind=engine)

    # Get valid chromosomes
    chroms = ParseUtils.get_chroms(session)

    # Get samples
    samples = _get_samples(session, samples_file)

    # Get source
    source = Source()
    source.name = source_name
    ParseUtils.upsert_source(session, source)
    source = ParseUtils.get_source(session, source_name)

    # This is ABSOLUTELY necessary to prevent MySQL from crashing!
    session.close()
    engine.dispose()

    # Parse metadata
    metadata = _parse_metadata(genome, metadata_file)

    exit(0)
    print(len(metadata))

    # Filter metadata (i.e. keep the best BED file per experiment)
    filtered_metadata = _filter_metadata(metadata)

    print(len(filtered_metadata))

    exit(0)

    # Group metadata by experiment and target
    grouped_metadata = _group_metadata(filtered_metadata)
    exit(0)

    # Parallelize inserts to the database
    _process_data_in_chunks(dummy_file, _insert_data_in_chunks, threads)

    # Dispose session
    Session.remove()

    # # Remove downloaded file
    # if os.path.exists(dummy_file):
    #     os.remove(dummy_file)

def _download_metadata(genome, feat_type, dummy_dir="/tmp/"):

    # Initialize
    url = "https://www.encodeproject.org/metadata/type=Experiment&status=released"

    # Fix hg38
    if genome == "hg38":
        genome = "GRCh38"

    # Add experiment to URL
    if feat_type == "atac-seq":
        experiment_url = "&assay_title=ATAC-seq&assay_slims=DNA%2Baccessibility&files.file_type=bam"
    elif feat_type == "accessibility":
        experiment_url = "&assay_title=DNase-seq&assay_title=FAIRE-seq&assay_slims=DNA%2Baccessibility&files.file_type=bed%2BnarrowPeak"
    elif feat_type == "histone":
        experiment_url = "&assay_title=Histone%2BChIP-seq&assay_slims=DNA%2Bbinding&files.file_type=bed%2BnarrowPeak"
    else:
        # TF
        experiment_url = "&assay_title=TF%2BChIP-seq&assay_slims=DNA%2Bbinding&files.file_type=bed%2BnarrowPeak"
    url += experiment_url

    # Add genome assembly
    url += "&assembly=%s" % genome

    # Download data
    metadata_file = os.path.join(url, "metadata.tsv")
    dummy_file = os.path.join(dummy_dir, "metadata.%s.%s.tsv" % (genome, feat_type))
    if not os.path.exists(dummy_file):
        f = urlretrieve(metadata_file, dummy_file)

    return(dummy_file)

def _get_samples(session, file_name):

    # Initialize
    samples = {}

    # For each line...
    for line in ParseUtils.parse_tsv_file(file_name):

        # If add...
        if line[3] == "Yes":

            # Get sample
            sample = Sample()
            sample.name = line[2]
            sample.treatment = False
            if line[4] == "Yes":
                sample.treatment = True
            sample.cell_line = False
            if line[5] == "Yes":
                sample.cell_line = True
            sample.cancer = False
            if line[6] == "Yes":
                sample.cancer = True

            # Upsert sample
            ParseUtils.upsert_sample(session, sample)

            # Get sample ID
            sample = ParseUtils.get_sample(session, sample.name, sample.X, sample.Y, sample.treatment, sample.cell_line, sample.cancer)

            # Add sample
            # samples.setdefault(line[0], sample.uid)
            samples.setdefault(line[0], None)

    return(samples)

def _parse_metadata(genome, metadata_file, test=True):

    # Initialize
    metadata_objects = set()
    accession_idx = None

    # For each line...
    for line in ParseUtils.parse_tsv_file(metadata_file):

        # If first line...
        if accession_idx is not None:

            # Initialize
            accession = line[accession_idx]
            biosample = line[biosample_idx]
            download_url = line[download_idx]
            experiment_accession = line[experiment_acc_idx]
            experiment_type = line[experiment_type_idx]
            experiment_target = line[experiment_target_idx]
            genome_assembly = line[genome_assembly_idx]
            output_format = line[output_format_idx]
            output_type = line[output_type_idx]
            status = line[status_idx]
            treatments = line[treatment_idx]

            if test:
                if biosample not in samples:
                    print(biosample)

            # ENCODE metadata object
            metadata_object = EncodeMetadata(accession, biosample, download_url, experiment_accession, experiment_type, experiment_target, genome_assembly, output_format, output_type, status, treatments)

            # Add metadata
            if metadata_object.genome_assembly == genome and metadata_object.status == "released" and type(metadata_object.treatments) is float:
                metadata_objects.add(metadata_object)

        else:
            # Get indices
            accession_idx = line.index("File accession")
             = line.index("Biosample term name")
            download_idx = line.index("File download URL")
            experiment_acc_idx = line.index("Experiment accession")
            experiment_type_idx = line.index("Assay")
            experiment_target_idx = line.index("Experiment target")
            genome_assembly_idx = line.index("Assembly")
            output_format_idx = line.index("File format")
            output_type_idx = line.index("Output type")
            status_idx = line.index("File Status")
            treatment_idx = line.index("Biosample treatments")

    return(metadata_objects)

def _filter_metadata(metadata):

    # Initialize
    grouped_metadata = {}
    filtered_metadata = set()
    output_types = ["optimal idr thresholded peaks", "pseudoreplicated idr thresholded peaks", "peaks", "alignments"]

    # For each ENCODE metadata object...
    for metadata_object in metadata:

        # Group metadata by experiment accession
        grouped_metadata.setdefault(metadata_object.experiment_accession, [])
        grouped_metadata[metadata_object.experiment_accession].append(metadata_object)

    # For each experiment accession...
    for experiment_acc in grouped_metadata:

        # Initialize
        exit_loop = False

        # For each output type...
        for out_type in output_types:

            # For each ENCODE metadata object...
            for metadata_object in grouped_metadata[experiment_acc]:

                if metadata_object.output_type == out_type:
                    filtered_metadata.add(metadata_object)
                    exit_loop = True

            if exit_loop:
                break

    return(filtered_metadata)

def _group_metadata(metadata):

    # Initialize
    grouped_metadata = {}

    # For each ENCODE metadata object...
    for metadata_object in metadata:

        # Group metadata
        grouped_metadata.setdefault()

    return(grouped_metadata)

        # # Skip tagged samples
        # if tag:
        #     print(": "\
        #         .join(
        #             [
        #                 os.path.basename(__file__),
        #                 "warning",
        #                 "use of protein tag",
        #                 "\"%s\" (\"%s\")" % (
        #                     accession,
        #                     tag
        #                 )
        #             ]
        #         )
        #     )
        # # Skip treated samples
        # if treatment:
        #     print(": "\
        #         .join(
        #             [
        #                 os.path.basename(__file__),
        #                 "warning",
        #                 "treated sample",
        #                 "\"%s\" (\"%s\")" % (
        #                     accession,
        #                     treatment
        #                 )
        #             ]
        #         )
        #     )
        #     continue
        # # This is a released sample!
        # if assembly == genome and status == "released":
        #     # Skip sample
        #     if not samples[biosample]["add"]: continue
        #     # Get metadata
        #     if os.path.exists(
        #         os.path.join(
        #             data_dir,
        #             "%s.bed.gz" % accession
        #         )
        #     ):
        #         k = (
        #             experiment_type,
        #             experiment_target
        #         )
        #         metadata.setdefault(k, [])
        #         metadata[k].append((accession, biosample))

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
        region.chrom = line[1]
        region.start = int(line[2])
        region.end = int(line[3])
        region.bin = assign_bin(region.start, region.end)

        # Ignore non-standard chroms, scaffolds, etc.
        if region.chrom not in chroms:
            continue

        # Upsert region
        _upsert_region(session, region)

        # Get region ID
        region = _get_region(session, region.chrom, region.start, region.end, region.strand)

        # Get conservation
        conservation = Conservation()
        conservation.regionID = region.uid
        conservation.sourceID = source.uid
        conservation.score = float(line[6])

        # Upsert conservation
        _upsert_conservation(session, conservation)

    # This is ABSOLUTELY necessary to prevent MySQL from crashing!
    session.close()
    engine.dispose()

#-------------#
# Main        #
#-------------#

if __name__ == "__main__": main()