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
from GUD import GUDglobals
from GUD.ORM.dna_accessibility import DNAAccessibility
from GUD.ORM.experiment import Experiment
from GUD.ORM.histone_modification import HistoneModification
from GUD.ORM.region import Region
from GUD.ORM.sample import Sample
from GUD.ORM.source import Source
from GUD.ORM.tf_binding import TFBinding
from . import _get_chroms, _get_db_name, _get_experiment, _get_region, _get_sample, _get_source, _initialize_gud_db, _initialize_engine_session, _process_data_in_chunks, _upsert_experiment, _upsert_region, _upsert_sample, _upsert_source, _upsert_tf

usage_msg = """
usage: %s --genome STR [-h] [options]
""" % os.path.basename(__file__)

help_msg = """%s
inserts features from ENCODE (Encyclopedia of DNA Elements)
into GUD.

  --genome STR        genome assembly
  --samples FILE      ENCODE samples (manually-curated)
  --feature STR       type of genomic feature

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
# Classes     #
#-------------#

class EncodeMetadata:

    def __init__(self, metadata=[]):
        """
        0 = File_accession
        1 = File_format
        2 = Output_type
        3 = Experiment_accession
        4 = Assay
        5 = Biosample_term_id
        6 = Biosample_term_name
        7 = Biosample_type
        8 = Biosample_organism
        9 = Biosample_treatments
        10 = Biosample_treatments_amount
        11 = Biosample_treatments_duration
        12 = Biosample_genetic_modifications_methods
        13 = Biosample_genetic_modifications_categories
        14 = Biosample_genetic_modifications_targets
        15 = Biosample_genetic_modifications_gene_targets
        16 = Biosample_genetic_modifications_site_coordinates
        17 = Biosample_genetic_modifications_zygositExperiment_target
        18 = Library_made_from
        19 = Library_depleted_in
        20 = Library_extraction_method
        21 = Library_lysis_method
        22 = Library_crosslinking_method
        23 = Library_strand_specific
        24 = Experiment_date_released
        25 = Project
        26 = Rbns_protein_concentration
        27 = Library_fragmentation_method
        28 = Library_size_range
        29 = Biological_replicates
        30 = Technical_replicate
        31 = Read_length
        32 = Mapped_read_length
        33 = Run_type
        34 = Paired_end
        35 = Paired_with
        36 = Derived_from
        37 = Size
        38 = Lab
        39 = md5sum
        40 = dbxrefs
        41 = File_download_URL
        42 = AssemblGenome_annotation
        43 = Platform
        44 = Controlled_by
        45 = File_Status
        46 = s3_uri
        47 = Audit_WARNING
        48 = Audit_INTERNAL_ACTION
        49 = Audit_NOT_COMPLIANT
        50 = Audit_ERROR
        """

        self.File_accession = metadata[0]
        self.File_format = metadata[1]
        self.Output_type = metadata[2]
        self.Experiment_accession = metadata[3]
        self.Assay = metadata[4]
        self.Biosample_term_id = metadata[5]
        self.Biosample_term_name = metadata[6]
        self.Biosample_type = metadata[7]
        self.Biosample_organism = metadata[8]
        self.Biosample_treatments = metadata[9]
        self.Biosample_treatments_amount = metadata[10]
        self.Biosample_treatments_duration = metadata[11]
        self.Biosample_genetic_modifications_methods = metadata[12]
        self.Biosample_genetic_modifications_categories = metadata[13]
        self.Biosample_genetic_modifications_targets = metadata[14]
        self.Biosample_genetic_modifications_gene_targets = metadata[15]
        self.Biosample_genetic_modifications_site_coordinates = metadata[16]
        self.Biosample_genetic_modifications_zygositExperiment_target = metadata[17]
        self.Library_made_from = metadata[18]
        self.Library_depleted_in = metadata[19]
        self.Library_extraction_method = metadata[20]
        self.Library_lysis_method = metadata[21]
        self.Library_crosslinking_method = metadata[22]
        self.Library_strand_specific = metadata[23]
        self.Experiment_date_released = metadata[24]
        self.Project = metadata[25]
        self.Rbns_protein_concentration = metadata[26]
        self.Library_fragmentation_method = metadata[27]
        self.Library_size_range = metadata[28]
        self.Biological_replicates = metadata[29]
        self.Technical_replicate = metadata[30]
        self.Read_length = metadata[31]
        self.Mapped_read_length = metadata[32]
        self.Run_type = metadata[33]
        self.Paired_end = metadata[34]
        self.Paired_with = metadata[35]
        self.Derived_from = metadata[36]
        self.Size = metadata[37]
        self.Lab = metadata[38]
        self.md5sum = metadata[39]
        self.dbxrefs = metadata[40]
        self.File_download_URL = metadata[41]
        self.AssemblGenome_annotation = metadata[42]
        self.Platform = metadata[43]
        self.Controlled_by = metadata[44]
        self.File_Status = metadata[45]
        self.s3_uri = metadata[46]
        self.Audit_WARNING = metadata[47]
        self.Audit_INTERNAL_ACTION = metadata[48]
        self.Audit_NOT_COMPLIANT = metadata[49]
        self.Audit_ERROR = metadata[50]

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
    if not args.genome or not args.samples or not args.feature:
        error = ["%s\n%s" % (usage_msg, os.path.basename(__file__)), "error", "argument \"--genome\" \"--samples\" \"--feature\" are required\n"]
        print(": ".join(error))
        exit(0)

    # Check for invalid feature
    valid_features = ["accessibility", "histone", "tf"]
    if args.feature not in valid_features:
        error = ["%s\n%s" % (usage_msg, os.path.basename(__file__)), "error", "argument \"--feature\"", "invalid choice", "\"%s\" (choose from" % args.feature, "%s)\n" % " ".join(["\"%s\"" % i for i in feats])]
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

    # Insert ENCODE data
    encode_to_gud(args.user, args.pwd, args.host, args.port, args.db, args.genome, args.samples, args.feature, args.merge, args.dummy_dir, args.threads)

def encode_to_gud(user, pwd, host, port, db, genome, samples_file, feat_type, merge=False, dummy_dir="/tmp/", threads=1):
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
    db_name = _get_db_name(user, pwd, host, port, db)

    # If database does not exist...
    if not database_exists(db_name):
        _initialize_gud_db(user, pwd, host, port, db, genome)

    # Get engine/session
    engine, Session = _initialize_engine_session(db_name)
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
    chroms = _get_chroms(session)

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

    # Parse metadata
    metadata = _parse_metadata()

    # Parallelize inserts to the database
    _process_data_in_chunks(dummy_file, _insert_data_in_chunks, threads)

    # Dispose session
    Session.remove()

    # # Remove downloaded file
    # if os.path.exists(dummy_file):
    #     os.remove(dummy_file)

def _download_metadata(genome, feat_type, dummy_dir="/tmp/"):

    # Initialize
    url = "https://www.encodeproject.org/metadata/type=Experiment"

    # Fix hg38
    if genome == "hg38":
        genome = "GRCh38"

    # Add experiment to URL
    if feat_type == "histone":
        experiment_url = "&assay_title=Histone%2BChIP-seq"
    url += experiment_url

    # Add organism to URL
    if genome == "hg19" or genome == "GRCh38":
        organism_url = "&replicates.library.biosample.donor.organism.scientific_name=Homo%2Bsapiens"
    url += organism_url

    # Add genome assembly
    assembly_url = "&assembly=%s" % genome
    url += assembly_url

    # Add extra to URL
    extra_url = "&files.file_type=bed%2BnarrowPeak&status=released"
    url += extra_url

    # Download data
    metadata_file = os.path.join(url, "metadata.tsv")
    dummy_file = os.path.join(dummy_dir, "metadata.tsv")
    if not os.path.exists(dummy_file):
        f = urlretrieve(metadata_file, dummy_file)

    return(dummy_file)

def _get_samples(session, file_name):

    # Initialize
    samples = {}

    # For each line...
    for line in GUDglobals.parse_tsv_file(file_name):

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
            _upsert_sample(session, sample)

            # Get sample ID
            sample = _get_sample(session, sample.name, sample.X, sample.Y, sample.treatment, sample.cell_line, sample.cancer)

            # Add sample
            samples.setdefault(line[0], sample.uid)

    return(samples)

def _parse_metadata(file_name):

    # Initialize
    metadata = {}

# For each line...
    for line in GUDglobals.parse_tsv_file(
        metadata_file
    ):
        # If first line...
        if accession_idx is None:
            accession_idx = line.index(
                "File accession"
            )
            assembly_idx = line.index(
                "Assembly"
            )
            biosample_idx = line.index(
                "Biosample term name"
            )
            experiment_type_idx = line.index(
                "Assay"
            )
            experiment_target_idx = line.index(
                "Experiment target"
            )
            treatment_idx = line.index(
                "Biosample treatments"
            )
            status_idx = line.index(
                "File Status"
            )
            continue
        # Initialize
        accession = \
            line[accession_idx]
        experiment_type = \
            line[experiment_type_idx]
        biosample = \
            line[biosample_idx]
        tag = None
        experiment_target = None
        m = re.search(
            "^(3xFLAG|eGFP)?-?(.+)-(human|mouse)$",
            line[experiment_target_idx]
        )
        if m:
            tag = m.group(1)
            experiment_target = m.group(2)
        treatment = \
            line[treatment_idx]
        assembly = line[assembly_idx]
        status = line[status_idx]
        # Skip tagged samples
        if tag:
            print(": "\
                .join(
                    [
                        os.path.basename(__file__),
                        "warning",
                        "use of protein tag",
                        "\"%s\" (\"%s\")" % (
                            accession,
                            tag
                        )
                    ]
                )
            )
        # Skip treated samples
        if treatment:
            print(": "\
                .join(
                    [
                        os.path.basename(__file__),
                        "warning",
                        "treated sample",
                        "\"%s\" (\"%s\")" % (
                            accession,
                            treatment
                        )
                    ]
                )
            )
            continue
        # This is a released sample!
        if assembly == genome and status == "released":
            # Skip sample
            if not samples[biosample]["add"]: continue
            # Get metadata
            if os.path.exists(
                os.path.join(
                    data_dir,
                    "%s.bed.gz" % accession
                )
            ):
                k = (
                    experiment_type,
                    experiment_target
                )
                metadata.setdefault(k, [])
                metadata[k].append((accession, biosample))


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