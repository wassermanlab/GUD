#!/usr/bin/env python

import argparse
from binning import assign_bin
import getpass
import os
from pybedtools import BedTool, cleanup, set_tempdir
import re
import shutil
import subprocess

# Import from GUD module
from GUD import GUDglobals
from GUD.ORM.sample import Sample
from GUD.ORM.tf_binding import TFBinding
#from .bed2gud import bed_to_gud_db
from .initialize import initialize_gud_db
from . import _get_chroms, _get_db_name, _get_experiment, _get_session, _get_source

usage_msg = """
usage: %s --data-dir DIR --samples FILE [-h] [options]
""" % os.path.basename(__file__)

help_msg = """%s
inserts ReMap transcription factor ChIP-seq data into GUD.

  --data-dir DIR      directory where data was downloaded
  --samples FILE      ENCODE samples (manually-curated)

optional arguments:
  -h, --help          show this help message and exit
  -c, --cluster       cluster genome regions by UCSC's regCluster 
                      (default = False)
  --dummy-dir DIR     dummy directory (default = "/tmp/")
  --source STR        source name (default = "ReMap")

mysql arguments:
  -d STR, --db STR    database name (default = "%s")
  -H STR, --host STR  host name (default = "localhost")
  -p STR, --pwd STR   password (default = ignore this option)
  -P STR, --port STR  port number (default = %s)
  -u STR, --user STR  user name (default = current user)
""" % (usage_msg, GUDglobals.db_name, GUDglobals.db_port)

#-------------#
# Functions   #
#-------------#

def parse_args():
    """
    This function parses arguments provided via the command line and returns an {argparse} object.
    """

    parser = argparse.ArgumentParser(add_help=False)

    # Mandatory args
    parser.add_argument("--data-dir")
    parser.add_argument("--samples")

    # Optional args
    optional_group = parser.add_argument_group("optional arguments")
    optional_group.add_argument("-h", "--help", action="store_true")
    optional_group.add_argument("-c", "--cluster", action="store_true")
    optional_group.add_argument("--dummy-dir", default="/tmp/")
    optional_group.add_argument("--source", default="ReMap")

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
    if not args.data_dir or not args.samples:
        error = ["%s\n%s" % (usage_msg, os.path.basename(__file__)), "error", "arguments \"--data-dir\" \"--samples\" are required\n"]
        print(": ".join(error))
        exit(0)

    # Check MySQL password
    if not args.pwd:
        args.pwd = ""

def main():

    # Parse arguments
    args = parse_args()

    # Insert ENCODE data
    remap_to_gud(args.user, args.pwd, args.host, args.port, args.db, args.data_dir, args.samples, args.cluster, args.dummy_dir, args.source)

def remap_to_gud(user, pwd, host, port, db, data_dir, samples_file, cluster=False, dummy_dir="/tmp/", source_name="ReMap"):

    # Initialize
    global session
    table = TFBinding()
    experiment_type = "ChIP-seq"
    set_tempdir(dummy_dir) # i.e. for pyBedTools

    # Get DB name
    db_name = _get_db_name(user, pwd, host, port, db)

    # Get session
    session = _get_session(db_name)

    # Get source
    global sou
    sou = _get_source(source_name)

    # Get experiment
    global exp
    exp = _get_experiment(session, experiment_type)

    # Get samples
    global samples
    samples = _get_samples(sample_file)

    # Get sample to uid
    global sample2uid
    sample2uid = _get_sample_to_uid(session, samples)

    if not engine.has_table(table.__tablename__):
        # Create table
        table.__table__.create(bind=engine)

    # Unwind ChIP-seq data from 1x BED per TF to 1x BED per TF per sample
    _unwind_bed_files(data_dir, dummy_dir)

    # For each directory...
    for d in os.listdir(dummy_dir):
        # Initialize
        dummy_file = os.path.join(dummy_dir, "%s.bed" % d)
        exp_dummy_dir = os.path.join(dummy_dir, d)
        # Skip if not directory
        if not os.path.isdir(exp_dummy_dir):
            continue
        # Skip if completed
        completed_file = os.path.join(dummy_dir, "%s.completed" % directory)
        if os.path.exists(completed_file):
            continue
        # Sort BED files
        _sort_bed_files(exp_dummy_dir)
        # Cluster BED files
        if cluster:
            _cluster(exp_dummy_dir)
        # Do not cluster
        else:
            _do_not_cluster(exp_dummy_dir)

        # Completed
        open(completed_file, "a").close()
#        # Do not cluster
#        else:
#            # For each accession, biosample...
#            for accession, biosample in metadata[k]:
#                # If BED file exists...
#                bed_file = os.path.join(
#                    exp_dummy_dir, "%s.bed" % accession)
#                if os.path.exists(bed_file):
#                    # Initialize
#                    histone_type = None
#                    tf_name = None
#                    if feat_type == "histone":
#                        histone_type = experiment_target
#                    if feat_type == "tf":
#                        tf_name = experiment_target
#                    # Insert BED file to GUD database
#                    bed_to_gud_db(
#                        user,
#                        pwd,
#                        host,
#                        port,
#                        db,
#                        bed_file,
#                        feat_type,
#                        exp.name,
#                        samples[biosample]["cell_or_tissue"],
#                        sou.name,
#                        histone_type,
#                        None,
#                        tf_name,
#                        samples[biosample]["cancer"],
#                        samples[biosample]["cell_line"],
#                        samples[biosample]["treatment"]
#                    )
#        # Remove dummy dir
#        if os.path.isdir(exp_dummy_dir): shutil.rmtree(exp_dummy_dir)

def _get_samples(file_name):

    # Initialize
    samples = {}

    # For each line...
    for line in GUDglobals.parse_tsv_file(file_name):
        # Initialize 
        cell_type = line[0]
        cell_or_tissue = line[6]
        if line[7] == "Yes": add = True
        else: add = False
        if line[8] == "Yes": treatment = True
        else: treatment = False
        if line[9] == "Yes": cell_line = True
        else: cell_line = False
        if line[10] == "Yes": cancer = True
        else: cancer = False
        # Add sample
        samples.setdefault(cell_type, {"cell_or_tissue": cell_or_tissue, "treatment": treatment, "cell_line": cell_line, "cancer": cancer, "add": add})

    return(samples)

def _get_sample_to_uid(session, samples):

    # Initialize
    sample2uid = {}

    for s in samples:
        # Initialize
        add = samples[s]["add"]
        cell_or_tissue = samples[s]["cell_or_tissue"]
        treatment = samples[s]["treatment"]
        cell_line = samples[s]["cell_line"]
        cancer = samples[s]["cancer"]
        # Skip sample
        if not add:
            continue
        # Get sample
        sample = Sample()
        if sample.is_unique(session, cell_or_tissue, treatment, cell_line, cancer):
            sample.name = cell_or_tissue
            sample.treatment = treatment
            sample.cell_line = cell_line
            sample.cancer = cancer
            session.add(sample)
            session.commit()
        sam = sample.select_unique(session, cell_or_tissue, treatment, cell_line, cancer)
        sample2uid.setdefault(s, sam.uid)

    return(sample2uid)

def _unwind_bed_files(data_dir, dummy_dir="/tmp/"):

    # For each file...
    for bed_file in os.listdir(data_dir):
        # Skip non-BED files
        if not bed_file.endswith(".bed") and not bed_file.endswith(".bed.gz"):
            continue
        # Get TF name
        m = re.search("^remap\d{4}_(.+)_all_.+.bed", bed_file)
        if m:
            # Initialize
            exp_dummy_dir = os.path.join(dummy_dir, "ChIP-seq.%s" % m.group(1))
            # Skip if already done
            if os.path.isdir(exp_dummy_dir):
                continue
            # Create dummy dir
            if not os.path.isdir(exp_dummy_dir):
                os.mkdir(exp_dummy_dir)
            # For each line...
            bed_file = os.path.join(data_dir, bed_file)
            for line in GUDglobals.parse_tsv_file(bed_file):
                # Initialize
                d, tf, s = line[3].split(".")
                tf = tf.split("_")
                s = s.split("_")
                # Fix sample
                if s[0] == "LNCAPABL":
                    s[0] == "LNCAP"
                # If sample exists...
                if s[0] in samples and tf[0] == m.group(1):
                    # Write
                    dummy_file = os.path.join(exp_dummy_dir, "%s.bed" % s[0])
                    GUDglobals.write(dummy_file, "\t".join(line))

def _sort_bed_files(directory):

    # For each BED file...
    for bed_file in os.listdir(directory):
        # Initialize
        lines = []
        # Load BED file
        bed = BedTool(os.path.join(directory, bed_file))
        # For each line...
        for line in bed:
            # If BED7+...
            if len(line) >= 7:
                extra = [str(int(l[6])+1), "0,0,0"]
                lines.append("\t".join(line[:7] + extra))
        # Load BED from string
        bed = BedTool("\n".join(lines), from_string=True)
        # Sort BED
        bed = bed.sort()
        # Save BED
        bed.saveas(dummy_file)
        # Copy dummy BED to original
        shutil.copy(dummy_file, file_name)
        # Remove dummy BED
        os.remove(dummy_file)
    # Clean PyBedTools files
    cleanup(remove_all=True)

def _cluster(directory):

    # Initialize
    file_names = []

    # For each BED file...
    for bed_file in os.listdir(exp_dummy_dir):
        # Skip non-BED files
        if not bed_file.endswith(".bed"):
            continue
        # Skip regCluster BED file
        if bed_file == "regCluster.bed":
            continue
        # Add BED file
        file_names.append(os.path.join(exp_dummy_dir, bed_file))
    # If multiple BED files...
    if len(file_names) > 1:
        _regCluster(file_names, directory)
    # ... Else...
    else:
        _merge(file_names[0], directory)

def _regCluster(file_names, directory):

    # Initialize
    label2sample = {}

    # Get TF name
    m = re.search("ChIP-seq.(.+)$", directory)
    experiment_target = m.group(1)
    # Run regCluster
    table_file, cluster_file, bed_file = _run_regCluster(file_names, directory)
    # For each line...
    for line in GUDglobals.parse_tsv_file(table_file):
        m = re.search("%s/(\w+).bed$" % str(experiment_target), line[0])
        if m:
            label2sample.setdefault(line[-1], m.group(1))
    # Insert regions
    _insert_regions_from_bed(bed_file)
    # Get regions
    regions = _get_regions_from_bed(bed_file)
    # Insert features
    _insert_features(regions, label2sample)

def _run_regCluster(file_names, directory):

    # If BED list file does not exist...
    bed_files = os.path.join(directory, "files.txt")
    if not os.path.exists(bed_files):
        # For each BED file...
        for bed_file in file_names:
            # Skip non-BED files
            if not bed_file.endswith(".bed"):
                continue
            # Skip regCluster BED file
            if bed_file == "regCluster.bed":
                continue
            # Add file to list
            GUDglobals.write(bed_files, bed_file)
    # If table of tables file does not exist...
    table_file = os.path.join(directory, "tableOfTables.txt")
    if not os.path.exists(table_file):
        # Make table of tables
        cmd = "regClusterMakeTableOfTables uw01 %s %s" (bed_files, table_file)
        process = subprocess.call(cmd, shell=True)
    # If cluster file does not exist 
    cluster_file = os.path.join(directory, "regCluster.cluster")
    bed_file = os.path.join(directory, "regCluster.bed")
    if not os.path.exists("%s.cluster" % cluster_file):
        # Cluster
        cmd = "regCluster %s %s %s" (table_file, cluster_file, bed_file)
        process = subprocess.call(cmd, shell=True)

    return(table_file, cluster_file, bed_file)



def _insert_features(regions, label2sample=None):
        # For each line...
    for line in GUDglobals.parse_tsv_file(
        "%s.cluster" % cluster_file
    ):
        # Get region
        reg_uid = regions[int(line[0]) - 1]
        # Skip invalid region
        if not reg_uid: continue
        # Get sample
        sam_uid =\
            accession2sample[
                label2accession[line[-1]]
            ]
        # Get TF feature
        feat = TFBinding()
        is_unique = feat.is_unique(
            session,
            reg_uid,
            sam_uid,
            exp.uid,
            sou.uid,
            experiment_target
        )
        # Insert feature to GUD
        if is_unique:
            feat.regionID = reg_uid
            feat.sampleID = sam_uid
            feat.experimentID =\
                exp.uid
            feat.sourceID = sou.uid
            feat.tf = experiment_target
            session.add(feat)
            session.commit()


def _merge(file_name):

        # Get sample
        m = re.search("^%s/(.+).bed$" % directory, bed_file)
        sam_uid = accession2sample[m.group(1)]
        # Load BED
        a = BedTool(bed_file)
        # For line in BED...
        for l in a.merge():
            # Get coordinates
            chrom = l[0]
            start = int(l[1])
            end = int(l[2])
            # If valid chromosome...
            if chrom in chroms:
                # Get region
                region = Region()
                if region.is_unique(
                    session,
                    chrom,
                    start,
                    end
                ):
                    # Insert region
                    region.bin =\
                        assign_bin(
                            start,
                            end
                        )
                    region.chrom = chrom
                    region.start = start
                    region.end = end
                    session.add(region)
                    session.commit()
                reg = region.select_unique(
                    session,
                    chrom,
                    start,
                    end
                )
                regions.append(reg.uid)
        # Clean PyBedTools files
        cleanup(remove_all=True)
        # For each region...
        for reg_uid in regions:
            # Get TF feature
            feat = TFBinding()
            is_unique = feat.is_unique(
                session,
                reg_uid,
                sam_uid,
                exp.uid,
                sou.uid,
                experiment_target
            )
            # Insert feature to GUD
            if is_unique:
                feat.regionID = reg_uid
                feat.sampleID = sam_uid
                feat.experimentID =\
                    exp.uid
                feat.sourceID = sou.uid
                feat.tf = experiment_target
                session.add(feat)
                session.commit()

#-------------#
# Main        #
#-------------#

if __name__ == "__main__": main()