#!/usr/bin/env python

import argparse
from binning import assign_bin
from functools import partial
import getpass
from multiprocessing import Pool, cpu_count
from numpy import isnan
import os
import pickle
import re
import requests
import shutil
import subprocess
import sys
import warnings
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
usage: %s --genome STR --samples FILE --feature STR [-h] [options]
""" % os.path.basename(__file__)

help_msg = """%s
inserts features from ENCODE (Encyclopedia of DNA Elements)
into GUD.

  --genome STR        genome assembly
  --samples FILE      ENCODE samples (manually-curated)
  --feature STR       type of genomic feature ("atac-seq",
                      "accessibility", "histone" or "tf")
  --sample-type STR   restrict to samples of speficied type
                      ("cell" or "tissue"; default = ignore)

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
# Classes     #
#-------------#

class ENCODE:

    def __init__(self, accession, biosample_name, biosample_type, download_url,
        experiment_accession, experiment_type, experiment_target, genome_assembly,
        output_format, output_type, status, treatments, genetic_modifications):

        # Fix hg38
        if genome_assembly == "GRCh38":
            genome_assembly = "hg38"

        self.accession = accession
        self.biosample_name = biosample_name
        self.biosample_type = biosample_type
        self.download_url = download_url
        self.experiment_accession = experiment_accession
        self.experiment_type = experiment_type
        self.experiment_target = experiment_target
        self.genome_assembly = genome_assembly
        self.output_format = output_format
        self.output_type = output_type
        self.status = status
        self.treatments = treatments
        self.genetic_modifications = genetic_modifications

        # To be initialized later
        self._biosample_sex = None
        self._biosample_summary = None
        self.cancer = False
        self.cell_line = False

    @property
    def sex(self):
        """
        sex of biosample
        """
        return(self._biosample_sex)

    @sex.setter
    def sex(self, value):
        if value == "female" or value == "male":
            self._biosample_sex = str(value)

    @property
    def summary(self):
        """
        summary of biosample
        """
        return(self._biosample_summary)

    @summary.setter
    def summary(self, value):
        self._biosample_summary = str(value)

    @property
    def X(self):
        """
        number of X chromosomes
        """
        if self.sex == "female":
            return(2)
        elif self.sex == "male":
            return(1)
        else:
            return(None)

    @property
    def Y(self):
        """
        number of Y chromosomes
        """
        if self.sex == "female":
            return(0)
        elif self.sex == "male":
            return(1)
        else:
            return(None)

    @property
    def treatment(self):
        """
        is biosamble treated?
        """
        return(self.treatments is not None)

    @property
    def genetic_modification(self):
        """
        has biosample been genetically modified?
        """
        return(self.genetic_modifications is not None)

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
    parser.add_argument("--sample-type")

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
    valid_features = ["atac-seq", "accessibility", "histone", "tf"]
    if args.feature not in valid_features:
        error = ["%s\n%s" % (usage_msg, os.path.basename(__file__)), "error", "argument \"--feature\"", "invalid choice", "\"%s\" (choose from" % args.feature, "%s)\n" % " ".join(["\"%s\"" % i for i in valid_features])]
        print(": ".join(error))
        exit(0)

    # Check for sample types
    valid_sample_types = ["cell", "tissue"]
    if args.sample_type is not None:
        if args.sample_type not in valid_sample_types:
            error = ["%s\n%s" % (usage_msg, os.path.basename(__file__)), "error", "argument \"--sample-type\"", "invalid choice", "\"%s\" (choose from" % args.sample_type, "%s or ignore this option)\n" % " ".join(["\"%s\"" % i for i in valid_sample_types])]
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

    # Insert ENCODE data
    encode_to_gud(args.genome, args.samples, args.feature, args.sample_type, args.dummy_dir, args.remove, args.test, args.threads)

def encode_to_gud(genome, samples_file, feat_type, sample_type=None, dummy_dir="/tmp/", remove=False, test=False, threads=1):
    """
    e.g. python -m GUD.parsers.encode2gud --genome hg38 --samples ./samples/ENCODE.tsv --feature accessibility
    """

    # Initialize
    global chroms
    global encodes
    global engine
    global experiment
    global samples
    global Feature
    global Session

    # Testing
    if test:
        global current_process
        from multiprocessing import current_process

    # Create dummy dir
    subdir = "%s.%s" % (genome, os.path.basename(__file__))
    dummy_dir = os.path.join(dummy_dir, subdir)
    if not os.path.isdir(dummy_dir):
        os.makedirs(dummy_dir)

    # Download metadata
    metadata_file = _download_metadata(genome, feat_type, dummy_dir)

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
    if feat_type == "accessibility" or feat_type == "atac-seq":
        Feature = DNAAccessibility
    elif feat_type == "histone":
        Feature = HistoneModification
    else:
        Feature = TFBinding
    ParseUtils.create_table(Feature)

    # Start a new session
    session = Session()

    # Get valid chromosomes
    chroms = ParseUtils.get_chroms(session)

    # Get samples
    samples = _get_samples(session, samples_file)

    # This is ABSOLUTELY necessary to prevent MySQL from crashing!
    session.close()
    engine.dispose()

    # Unless pickle file exists
    pickle_file = "%s.pickle" % metadata_file
    if not os.path.exists(pickle_file):

        # Parse metadata
        # Add biosample summary
        encodes = _add_biosample_summary(_parse_metadata(genome, metadata_file))
        handle = ParseUtils._get_file_handle(pickle_file, "wb")
        pickle.dump(encodes, handle)

    else:

        handle = ParseUtils._get_file_handle(pickle_file, "rb")
        encodes = pickle.load(handle)

    # Add biosample info (i.e. sex, cancer, cell line)
    encodes = _add_biosample_info(encodes)

    # Filter ENCODE accessions (i.e. for each experiment, keep the accession from the best type of file)
    # Group ENCODE accessions by experiment target and type
    grouped_accessions = _group_ENCODE_accessions(_filter_ENCODE_accessions(feat_type, sample_type))

    # For each experiment target/type...
    for experiment_target, experiment_type in sorted(grouped_accessions):

        # Beware, for this should not be possible!
        if experiment_target is not None:
            if feat_type != "histone" and feat_type != "tf":
                continue
        else:
            if feat_type == "histone" or feat_type == "tf":
                continue

        # Start a new session
        session = Session()

        # Get experiment
        experiment = Experiment()
        experiment.name = experiment_type
        ParseUtils.upsert_experiment(session, experiment)
        experiment = ParseUtils.get_experiment(session, experiment_type)

        # This is ABSOLUTELY necessary to prevent MySQL from crashing!
        session.close()
        engine.dispose()

        # Prepare data
        subgrouped_accessions = grouped_accessions[(experiment_target, experiment_type)]
        data_file = _preprocess_data(subgrouped_accessions, dummy_dir, test, threads)

        # Split data
        data_files = _split_data(data_file, threads)
        
        # Insert samples/sources
        _insert_samples_and_sources(subgrouped_accessions)

        # Parallelize inserts to the database
        ParseUtils.insert_data_files_in_parallel(data_files, partial(_insert_data_file, test=test), threads)

    # Remove files
    if remove:
        shutil.rmtree(dummy_dir)

    # Dispose session
    Session.remove()

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
        urlretrieve(metadata_file, dummy_file)

    return(dummy_file)

def _get_samples(session, file_name):

    # Initialize
    samples = {}

    # For each line...
    for line in ParseUtils.parse_tsv_file(file_name):

        # Initialize
        sample_name = line[0]
        cell_line = False
        cancer = False
        sex = None
        add = False

        if line[2] == "Yes":
            cell_line = True
        if line[3] == "Yes":
            cancer = True
        if line[4] == "female" or line[4] == "male":
            sex = line[4]
        if line[5] == "Yes":
            add = True

        # Add sample
        samples.setdefault(sample_name, (cell_line, cancer, sex, add))

    return(samples)

def _parse_metadata(genome, metadata_file):

    # Initialize
    accession_idx = None
    encode_objects = {}
    regexp = re.compile("^(3xFLAG|eGFP)?-?(.+)-(human|mouse|dmelanogaster)$")

    # For each line...
    for line in ParseUtils.parse_tsv_file(metadata_file):

        # If first line...
        if accession_idx is not None:

            # Initialize
            accession = line[accession_idx]
            biosample_name = line[biosample_name_idx]
            biosample_type = line[biosample_type_idx]
            download_url = line[download_idx]
            experiment_accession = line[experiment_acc_idx]
            experiment_type = line[experiment_type_idx]
            experiment_target = line[experiment_target_idx]
            if type(experiment_target) is float and isnan(experiment_target):
                experiment_target = None
            else:
                m = regexp.search(experiment_target)
                experiment_target = m.group(2)
            genome_assembly = line[genome_assembly_idx]
            output_format = line[output_format_idx]
            output_type = line[output_type_idx]
            status = line[status_idx]
            treatments = line[treatment_idx]
            if type(treatments) is float and isnan(treatments):
                treatments = None
            genetic_modifications = line[genetic_modifications_idx]
            if type(genetic_modifications) is float and isnan(genetic_modifications):
                genetic_modifications = None

            # Add ENCODE object
            encode = ENCODE(accession, biosample_name, biosample_type, download_url,
                experiment_accession, experiment_type, experiment_target, genome_assembly,
                output_format, output_type, status, treatments, genetic_modifications)
            if encode.genome_assembly == genome and encode.status == "released" and \
               not encode.treatment and not encode.genetic_modification:
                encode_objects.setdefault(encode.accession, encode)

        else:
            # Get indices
            accession_idx = line.index("File accession")
            biosample_name_idx = line.index("Biosample term name")
            biosample_type_idx = line.index("Biosample type")
            download_idx = line.index("File download URL")
            experiment_acc_idx = line.index("Experiment accession")
            experiment_type_idx = line.index("Assay")
            experiment_target_idx = line.index("Experiment target")
            genome_assembly_idx = line.index("File assembly")
            output_format_idx = line.index("File format")
            output_type_idx = line.index("Output type")
            status_idx = line.index("File Status") - 1
            treatment_idx = line.index("Biosample treatments")
            genetic_modifications_idx = line.index("Biosample genetic modifications methods")

    return(encode_objects)

def _add_biosample_summary(encode_objects):
    """
    https://www.encodeproject.org/help/rest-api/
    """

    # Initialize
    experiments = set()
    biosample_summaries = {}
    updated_encode_objects = {}
    headers = {"accept": "application/json"}

    # For each accession...
    for accession in encode_objects:

        # Get experiment
        encode = encode_objects[accession]
        experiments.add(encode.experiment_accession)

    # For each experiment...
    for experiment in experiments:

        # Get biosample summary
        url = "https://www.encodeproject.org/experiment/%s/?frame=object" % experiment
        response = requests.get(url, headers=headers)
        biosample_summary = response.json()["biosample_summary"]
        biosample_summaries.setdefault(experiment, biosample_summary)

    # For each accession...
    for accession in encode_objects:

        # Update ENCODE object
        encode = encode_objects[accession]
        encode.summary = biosample_summaries[encode.experiment_accession]
        updated_encode_objects.setdefault(accession, encode)

    return(updated_encode_objects)

def _add_biosample_info(encode_objects):

    # Initialize
    i_have_been_warned = False
    updated_encode_objects = {}

    # For each accession...
    for accession in encode_objects:

        encode = encode_objects[accession]

        # Warn me!
        if encode.summary not in samples:
            warnings.warn("missing sample: %s" % encode.summary, Warning, stacklevel=2)
            i_have_been_warned = True

        # Update ENCODE object
        elif samples[encode.summary][3]:
            encode.cancer = samples[encode.summary][1]
            encode.cell_line = samples[encode.summary][0]
            encode.sex = samples[encode.summary][2]
            updated_encode_objects.setdefault(accession, encode)

    if i_have_been_warned:
        error = ["%s" % os.path.basename(__file__), "error", "missing samples"]
        print(": ".join(error))
        exit(0)

    return(updated_encode_objects)

def _filter_ENCODE_accessions(feat_type, sample_type=None):

    # Initialize
    done = set()
    grouped_accessions = {}
    filtered_accessions = set()

    # Possible output files (from ENCODE pipelines)
    # https://www.encodeproject.org/atac-seq/
    # https://www.encodeproject.org/data-standards/dnase-seq/
    # https://www.encodeproject.org/chip-seq/histone/
    # https://www.encodeproject.org/chip-seq/transcription_factor/
    output_types = {
        "accessibility": ["peaks"],
        "atac-seq": [],
        "histone": ["peaks"],
        "tf": ["optimal idr thresholded peaks", "conservative idr thresholded peaks", "pseudoreplicated idr thresholded peaks", "peaks"]
    }

    # For each accession...
    for accession in encodes:

        # Group ENCODE objects by experiment accession
        encode = encodes[accession]
        if sample_type == "tissue" and encode.biosample_type != "tissue":
            continue
        elif sample_type == "cell" and encode.biosample_type == "tissue":
            continue
        grouped_accessions.setdefault(encode.experiment_accession, [])
        grouped_accessions[encode.experiment_accession].append(accession)

    # For each output type...
    for out_type in output_types[feat_type]:

        # For each experiment accession...
        for experiment_accession in grouped_accessions:        

            # Skip experiment
            if experiment_accession in done:
                continue

            # For each accession...
            for accession in grouped_accessions[experiment_accession]:

                if encodes[accession].output_type.lower() == out_type:
                    filtered_accessions.add(accession)
                    done.add(experiment_accession)

    return(filtered_accessions)

def _group_ENCODE_accessions(accessions):

    # Initialize
    grouped_accessions = {}

    # For each accession...
    for accession in accessions:

        # Initialize
        experiment_target = encodes[accession].experiment_target
        experiment_type = encodes[accession].experiment_type

        # Group metadata
        grouped_accessions.setdefault((experiment_target, experiment_type), set())
        grouped_accessions[(experiment_target, experiment_type)].add(accession)

    return(grouped_accessions)

def _preprocess_data(accessions, dummy_dir="/tmp/", test=False, threads=1):

    # Initialize
    dummy_files = []
    encode_objects = set()

    # Get label
    accession = next(iter(accessions))
    encode = encodes[accession]
    label = encode.experiment_type
    if encode.experiment_target is not None:
        label += ".%s" % encode.experiment_target

    # Skip if BED file exists
    bed_file = os.path.join(dummy_dir, "%s.bed" % label)
    if not os.path.exists(bed_file):

        # Skip if BED file exists
        dummy_file = os.path.join(dummy_dir, "dummy.bed")
        if not os.path.exists(dummy_file):

            # For each accession...
            for accession in accessions:
                encode_objects.add(encodes[accession])

            # Get ENCODE BED files
            pool = Pool(processes=threads)
            for download_file in pool.imap(partial(_download_ENCODE_bed_file, dummy_dir=dummy_dir, test=test), encode_objects):

                # Initialize
                m = re.search("\/(\w+).(bam|bed.gz)$", download_file)
                accession = m.group(1)

                # Concatenate
                cmd = "zless %s | awk '{print $1\"\t\"$2\"\t\"$3\"\t%s\t\"$7\"\t\"$10}' >> %s" % (download_file, accession, dummy_file)
                subprocess.call(cmd, shell=True)

                # Remove downloaded file
                os.remove(download_file)

            pool.close()

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

def _download_ENCODE_bed_file(encode, dummy_dir="/tmp/", test=False):

    # Initialize
    download_file = os.path.join(dummy_dir, encode.accession)

    # Testing
    if test:
        print(current_process().name)

    # DO NOT REMOVE THIS!
    # # Preprocess BAM data
    # if encode.output_format == "bam":

    #     # Skip if peaks file already exists
    #     peaks_file = os.path.join(dummy_dir, "%s_peaks.narrowPeak" % encode.accession)
    #     if not os.path.exists(peaks_file):

    #         # Download BAM file
    #         download_file += ".bam"
    #         if not os.path.exists(download_file):
    #            urlretrieve(encode.download_url, download_file)        

    #         # From "An atlas of chromatin accessibility in the adult human brain" (Fullard et al. 2018):
    #         # Peaks for OCRs were called using the model-based Analysis of ChIP-seq (MACS) (Zhang et al. 2008)
    #         # v2.1 (https://github.com/taoliu/MACS/). It models the shift size of tags and models local biases
    #         # in sequencability and mapability through a dynamic Poisson background model. We used the following
    #         # parameters (Kaufman et al. 2016): "--keep-dup all", "--shift -100", "--extsize 200", "--nomodel".
    #         cmd = "macs2 callpeak -t %s --keep-dup all --shift -100 --extsize 200 --nomodel --outdir %s -n %s" % (download_file, dummy_dir, encode.accession)
    #         process = subprocess.Popen([cmd], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    #         stdout, stderr = process.communicate()

    #         # Remove files
    #         os.remove(os.path.join(dummy_dir, "%s_peaks.xls" % encode.accession))
    #         os.remove(os.path.join(dummy_dir, "%s_summits.bed" % encode.accession))

    #         # Set peaks file as download file
    #         shutil.move(peaks_file, download_file)

    # else:

    # Download BED file
    download_file += ".bed.gz"
    if not os.path.exists(download_file):
        urlretrieve(encode.download_url, download_file)

    return(download_file)

def _insert_samples_and_sources(accessions):

    # Initialize
    session = Session()

    # For each accession...
    for accession in accessions:

        # Upsert sample
        sample = Sample()
        if not encodes[accession].summary:
            sample.name = encodes[accession].biosample_name
        else:
            sample.name = encodes[accession].summary
        sample.treatment = encodes[accession].treatment
        sample.cell_line = encodes[accession].cell_line
        sample.cancer = encodes[accession].cancer
        if encodes[accession].sex is not None:
            sample.X = encodes[accession].X
            sample.Y = encodes[accession].Y
        ParseUtils.upsert_sample(session, sample)

        # Upsert source
        source = Source()
        source.name = "ENCODE"
        # This is to satisfy ENCODE's citing guidelines:
        # https://www.encodeproject.org/help/citing-encode/
        source.source_metadata = "%s," % accession
        source.metadata_descriptor = "accession,"
        source.url = encodes[accession].download_url
        ParseUtils.upsert_source(session, source)

    # This is ABSOLUTELY necessary to prevent MySQL from crashing!
    session.close()
    engine.dispose()

def _insert_data_file(data_file, test=False):

    # Initialize
    session = Session()
    accession2sample = {}
    accession2source = {}

    # Testing
    if test:
        lines = 0
        print(current_process().name)

    # For each line...
    for line in ParseUtils.parse_tsv_file(data_file):

        # Initialize
        accession = line[3]

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
        if accession not in accession2sample:
            sample = Sample()
            if not encodes[accession].summary:
                sample.name = encodes[accession].biosample_name
            else:
                sample.name = encodes[accession].summary
            sample.treatment = encodes[accession].treatment
            sample.cell_line = encodes[accession].cell_line
            sample.cancer = encodes[accession].cancer
            if encodes[accession].sex is not None:
                sample.X = encodes[accession].X
                sample.Y = encodes[accession].Y
            sample = ParseUtils.get_sample(session, sample.name, sample.X, sample.Y, sample.treatment, sample.cell_line, sample.cancer)
            accession2sample.setdefault(accession, sample)
        else:
            sample = accession2sample[accession]

        # Get source
        if accession not in accession2source:
            source = Source()
            source.name = "ENCODE"
            source.source_metadata = "%s," % accession
            source.metadata_descriptor = "accession,"
            source.url = encodes[accession].download_url
            source = ParseUtils.get_source(session, source.name, source.source_metadata, source.metadata_descriptor, source.url)
            accession2source.setdefault(accession, source)
        else:
            source = accession2source[accession]

        # Upsert feature
        feature = Feature()
        feature.region_id = region.uid
        feature.sample_id = sample.uid
        feature.experiment_id = experiment.uid
        feature.source_id = source.uid
        feature.score = float(line[4])
        feature.peak = int(line[5])
        if Feature.__tablename__ == "dna_accessibility":
            ParseUtils.upsert_accessibility(session, feature)
        else:
            if Feature.__tablename__ == "histone_modifications":
                feature.histone_type = encodes[accession].experiment_target
                ParseUtils.upsert_histone(session, feature)
            else:
                feature.tf = encodes[accession].experiment_target
                ParseUtils.upsert_tf(session, feature)

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
