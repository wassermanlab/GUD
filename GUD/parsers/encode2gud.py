#!/usr/bin/env python

import argparse
from binning import assign_bin
from functools import partial
import getpass
from multiprocessing import Pool, cpu_count
from numpy import isnan
import os
from pybedtools import BedTool, cleanup, set_tempdir
import re
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
  --feature STR       type of genomic feature ("atac-seq"
                      "accessibility", "histone" or "tf")

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
    encode_to_gud(args.genome, args.samples, args.feature, args.dummy_dir, args.merge, args.test, args.threads)

def encode_to_gud(genome, samples_file, feat_type, dummy_dir="/tmp/", merge=False, test=False, threads=1):
    """
    python -m GUD.parsers.encode2gud --genome hg19 --samples --dummy-dir ./tmp/
    """

    # Initialize
    global chroms
    global engine
    global experiment
    global Feature
    global metadata
    global Session
    global samples
    global source
    source_name = "ENCODE"
    set_tempdir(dummy_dir) # i.e. for pyBedTools

    # Testing
    if test:
        global current_process
        from multiprocessing import current_process

    # Download metadata
    metadata_file = _download_metadata(genome, feat_type, dummy_dir)

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

    # Get source
    source = Source()
    source.name = source_name
    ParseUtils.upsert_source(session, source)
    sources = ParseUtils.get_source(session, source_name)
    source = next(iter(sources))

    # This is ABSOLUTELY necessary to prevent MySQL from crashing!
    session.close()
    engine.dispose()

    # Parse metadata
    metadata = _parse_metadata(genome, metadata_file)

    # Filter metadata (i.e. keep the best BED file per experiment)
    # Group metadata by experiment target and experiment type
    grouped_metadata = _group_metadata(_filter_metadata())

    # For each experiment target/type...
    for experiment_target, experiment_type in grouped_metadata:

        # Beware, for this is not possible!
        if experiment_target is not None:
            if feat_type != "histone" and feat_type != "tf":
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
        data_file = _preprocess_data(grouped_metadata[(experiment_target, experiment_type)], dummy_dir, merge, test, threads)

        # Split data
        data_files = _split_data(data_file, threads)

        # Parallelize inserts to the database
        ParseUtils.insert_data_files_in_parallel(data_files, partial(_insert_data_file, test=test), threads)

        # Remove data files
        if not test:
            if os.path.exists(data_file):
                os.remove(data_file)
            for data_file in data_files:
                if os.path.exists(data_file):
                    os.remove(data_file)

    # Remove downloaded file
    if os.path.exists(metadata_file) and not test:
        os.remove(metadata_file)

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
            samples.setdefault(line[0], sample.uid)

    return(samples)

def _parse_metadata(genome, metadata_file):

    # Initialize
    i_have_been_warned = False
    accession_idx = None
    metadata_objects = {}

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
            if type(experiment_target) is float and isnan(experiment_target):
                experiment_target = None
            genome_assembly = line[genome_assembly_idx]
            output_format = line[output_format_idx]
            output_type = line[output_type_idx]
            status = line[status_idx]
            treatments = line[treatment_idx]
            if type(treatments) is float and isnan(treatments):
                treatments = None

            # Warn me!
            if biosample not in samples:
                i_have_been_warned = True
                warnings.warn("missing sample: %s" % biosample, Warning, stacklevel=2)

            # ENCODE metadata object
            metadata_object = EncodeMetadata(accession, biosample, download_url, experiment_accession, experiment_type, experiment_target, genome_assembly, output_format, output_type, status, treatments)

            # Add metadata
            if metadata_object.genome_assembly == genome and metadata_object.status == "released" and not metadata_object.treatments:
                metadata_objects.setdefault(metadata_object.accession, metadata_object)

        else:
            # Get indices
            accession_idx = line.index("File accession")
            biosample_idx = line.index("Biosample term name")
            download_idx = line.index("File download URL")
            experiment_acc_idx = line.index("Experiment accession")
            experiment_type_idx = line.index("Assay")
            experiment_target_idx = line.index("Experiment target")
            genome_assembly_idx = line.index("Assembly")
            output_format_idx = line.index("File format")
            output_type_idx = line.index("Output type")
            status_idx = line.index("File Status")
            treatment_idx = line.index("Biosample treatments")

    if i_have_been_warned:
        error = ["%s" % os.path.basename(__file__), "error", "missing samples"]
        print(": ".join(error))
        exit(0)

    return(metadata_objects)

def _filter_metadata():

    # Initialize
    grouped_metadata = {}
    filtered_metadata = set()
    output_types = ["optimal idr thresholded peaks", "pseudoreplicated idr thresholded peaks", "peaks", "alignments"]

    # For each accession...
    for accession in metadata:

        # Group metadata by experiment accession
        metadata_object = metadata[accession]
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

        # Initialize
        experiment_target = metadata_object.experiment_target
        experiment_type = metadata_object.experiment_type

        # Group metadata
        grouped_metadata.setdefault((experiment_target, experiment_type), set())
        grouped_metadata[(experiment_target, experiment_type)].add(metadata_object)

    return(grouped_metadata)

def _preprocess_data(metadata_objects, dummy_dir="/tmp/", merge=False, test=False, threads=1):

    # Initialize
    dummy_files = []

    # Get label
    metadata_object = next(iter(metadata_objects))
    label = metadata_object.experiment_type
    if Feature.__tablename__ == "histone_modifications" or Feature.__tablename__ == "tf_binding":
        m = re.search("^(3xFLAG|eGFP)?-?(.+)-(human|mouse)$", metadata_object.experiment_target)
        label += ".%s" % m.group(2)

    # Skip if BED file exists
    bed_file = os.path.join(dummy_dir, "%s.bed" % label)
    if not os.path.exists(bed_file):

        # Skip if BED file exists
        dummy_file = os.path.join(dummy_dir, "dummy.bed")
        if not os.path.exists(dummy_file):

            # Get ENCODE BED files
            pool = Pool(processes=threads)
            for download_file in pool.imap(partial(_download_ENCODE_bed_file, dummy_dir=dummy_dir, test=test), metadata_objects):

                # Initialize
                m = re.search("\/(\w+).(bam|bed.gz)$", download_file)
                accession = m.group(1)

                # Concatenate
                cmd = "zless %s | awk '{print $1\"\t\"$2\"\t\"$3\"\t%s\"}' >> %s" % (download_file, accession, dummy_file)
                subprocess.call(cmd, shell=True)

                # Remove downloaded file
                os.remove(download_file)

            pool.close()

        # Add dummy file
        dummy_files.append(dummy_file)

        # Sort BED
        dummy_file = os.path.join(dummy_dir, "dummy.sorted.bed")
        if not os.path.exists(dummy_file):

            # UNIX sort is more efficient than bedtools
            cmd = "sort -T %s -k1,1 -k2,2n %s > %s" % (dummy_dir, dummy_files[0], dummy_file)
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
        a = BedTool(dummy_files[1])
        b = BedTool(dummy_files[-1])
        a.intersect(b, wa=True, wb=True, stream=True).saveas(bed_file)

        # Clean PyBedTools files
        cleanup(remove_all=True)

        # Remove ALL dummy files
        for dummy_file in dummy_files:
            os.remove(dummy_file)

    return(bed_file)

def _split_data(data_file, threads=1):

    # import math

    # Initialize
    split_files = []
    # data_file_dir = os.path.dirname(os.path.realpath(data_file))

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

    # # Get number of lines
    # process = subprocess.check_output(["wc -l %s" % data_file], shell=True)
    # m = re.search("(\d+)", str(process))
    # l = math.ceil(int(m.group(1))/(threads*10))

    # # Split
    # cmd = "split -l %s %s %s." % (l, data_file, data_file)
    # subprocess.call(cmd, shell=True)

    # # For each split file...
    # for split_file in os.listdir(data_file_dir):

    #     # Append split file
    #     split_file = os.path.join(data_file_dir, split_file)
    #     if os.path.abspath(data_file) in split_file:
    #         split_files.append(split_file)

    return(split_files)

def _download_ENCODE_bed_file(metadata_object, dummy_dir="/tmp/", test=False):

    # Initialize
    download_file = os.path.join(dummy_dir, metadata_object.accession)

    # Testing
    if test:
        print(current_process().name)

    # Preprocess BAM data
    if metadata_object.output_format == "bam":

        # Skip if peaks file already exists
        peaks_file = os.path.join(dummy_dir, "%s_peaks.narrowPeak" % metadata_object.accession)
        if not os.path.exists(peaks_file):

            # Download BAM file
            download_file += ".bam"
            if not os.path.exists(download_file):
               urlretrieve(metadata_object.download_url, download_file)        

            # From "An atlas of chromatin accessibility in the adult human brain" (Fullard et al. 2018):
            # Peaks for OCRs were called using the model-based Analysis of ChIP-seq (MACS) (Zhang et al. 2008)
            # v2.1 (https://github.com/taoliu/MACS/). It models the shift size of tags and models local biases
            # in sequencability and mapability through a dynamic Poisson background model. We used the following
            # parameters (Kaufman et al. 2016): "--keep-dup all", "--shift -100", "--extsize 200", "--nomodel".
            cmd = "macs2 callpeak -t %s --keep-dup all --shift -100 --extsize 200 --nomodel --outdir %s -n %s" % (download_file, dummy_dir, metadata_object.accession)
            process = subprocess.Popen([cmd], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            stdout, stderr = process.communicate()

            # Remove files
            os.remove(os.path.join(dummy_dir, "%s_peaks.xls" % metadata_object.accession))
            os.remove(os.path.join(dummy_dir, "%s_summits.bed" % metadata_object.accession))

            # Set peaks file as download file
            shutil.move(peaks_file, download_file)

    else:

        # Download BED file
        download_file += ".bed.gz"
        if not os.path.exists(download_file):
            urlretrieve(metadata_object.download_url, download_file)

    return(download_file)

def _insert_data_file(data_file, test=False):

    # Initialize
    session = Session()

    # Testing
    if test:
        lines = 0
        print(current_process().name)

    # For each line...
    for line in ParseUtils.parse_tsv_file(data_file):

        # Get region
        region = Region()
        region.chrom = line[-3]
        region.start = int(line[-2])
        region.end = int(line[-1])
        region.bin = assign_bin(region.start, region.end)

        # Ignore non-standard chroms, scaffolds, etc.
        if region.chrom not in chroms:
            continue

        # Get sample
        sample_name = metadata[line[3]].biosample

        # Ignore samples
        if sample_name not in samples:
            continue

        # Upsert region
        ParseUtils.upsert_region(session, region)

        # Get region ID
        region = ParseUtils.get_region(session, region.chrom, region.start, region.end, region.strand)

        # Get feature
        feature = Feature()
        feature.region_id = region.uid
        feature.sample_id = samples[sample_name]
        feature.experiment_id = experiment.uid
        feature.source_id = source.uid
        # Upsert accessibility
        if Feature.__tablename__ == "dna_accessibility":
            ParseUtils.upsert_accessibility(session, feature)
        else:
            # Get experiment target
            m = re.search("^(3xFLAG|eGFP)?-?(.+)-(human|mouse)$", metadata[line[3]].experiment_target)
            experiment_target = m.group(2)
            # Upsert histone
            if Feature.__tablename__ == "histone_modifications":
                feature.histone_type = experiment_target
                ParseUtils.upsert_histone(session, feature)
            # Upsert tf
            else:
                feature.tf = experiment_target
                ParseUtils.upsert_tf(session, feature)

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