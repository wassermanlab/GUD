#!/usr/bin/env python

import argparse
from binning import assign_bin
import getpass
import os
from pybedtools import BedTool, cleanup, set_tempdir
import re
import shutil
from sqlalchemy import create_engine
from sqlalchemy.orm import (
    scoped_session,
    sessionmaker
)
from sqlalchemy_utils import database_exists
import subprocess

# Import from GUD module
from GUD import GUDglobals
from GUD.ORM.experiment import Experiment
from GUD.ORM.region import Region
from GUD.ORM.sample import Sample
from GUD.ORM.source import Source
from GUD.ORM.tf_binding import TFBinding
#from .bed2gud import bed_to_gud_db
from .initialize import initialize_gud_db

usage_msg = """
usage: %s --data-dir DIR --samples FILE
                    [-h] [-c] [--dummy-dir DIR] [--source STR]
                    [-d STR] [-H STR] [-p STR] [-P STR] [-u STR]
""" % \
os.path.basename(__file__)

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
""" % \
(
    usage_msg,
    GUDglobals.db_name,
    GUDglobals.db_port
)


#-------------#
# Functions   #
#-------------#

def parse_args():
    """
    This function parses arguments provided via
    the command line and returns an {argparse}
    object.
    """

    parser = argparse.ArgumentParser(
        add_help=False,
    )

    # Mandatory args
    parser.add_argument("--data-dir")
    parser.add_argument("--samples")

    # Optional args
    optional_group = parser.add_argument_group(
        "optional arguments"
    )
    optional_group.add_argument(
        "-h", "--help",
        action="store_true"
    )
    optional_group.add_argument(
        "-c", "--cluster",
        action="store_true"
    )
    optional_group.add_argument(
        "--dummy-dir",
        default="/tmp/"
    )
    optional_group.add_argument(
        "--source",
        default="ReMap"
    )

    # MySQL args
    mysql_group = parser.add_argument_group(
        "mysql arguments"
    )
    mysql_group.add_argument(
        "-d", "--db",
        default=GUDglobals.db_name,
    )
    mysql_group.add_argument(
        "-H", "--host",
        default="localhost"
    )
    mysql_group.add_argument("-p", "--pwd")
    mysql_group.add_argument(
        "-P", "--port",
        default=GUDglobals.db_port
    )
    mysql_group.add_argument(
        "-u", "--user",
        default=getpass.getuser()
    )

    args = parser.parse_args()

    check_args(args)

    return args

def check_args(args):
    """
    This function checks an {argparse} object.
    """

    # Print help
    if args.help:
        print(help_msg)
        exit(0)

    # Check mandatory arguments
    if (
        not args.data_dir or \
        not args.samples
    ):
        print(": "\
            .join(
                [
                    "%s\n%s" % \
                        (
                            usage_msg,
                            os.path.basename(__file__)
                        ),
                    "error",
                    "arguments \"--data-dir\" \"--samples\" are required\n"
                ]
            )
        )
        exit(0)

    # Check MySQL password
    if not args.pwd:
        args.pwd = ""

def main():

    # Parse arguments
    args = parse_args()

    # Insert ENCODE data
    remap_to_gud(
        args.user,
        args.pwd,
        args.host,
        args.port,
        args.db,
        args.data_dir,
        args.samples,
        args.cluster,
        args.dummy_dir,
        args.source
    )

def remap_to_gud(user, pwd, host, port, db,
    data_dir, samples_file, cluster=False,
    dummy_dir="/tmp/", source_name="ReMap"):

    # Initialize
    samples = {}
    accessions2samples = {}
    table = TFBinding()
    experiment_type = "ChIP-seq"
    set_tempdir(dummy_dir)
    db_name = "mysql://{}:{}@{}:{}/{}".format(
        user, pwd, host, port, db
    )
    if not database_exists(db_name):
        initialize_gud_db(
            user,
            pwd,
            host,
            port,
            db,
            genome
        )
    session = scoped_session(sessionmaker())
    engine = create_engine(
        db_name,
        echo=False,
        pool_pre_ping=True
    )
    session.remove()
    session.configure(
        bind=engine,
        autoflush=False,
        expire_on_commit=False
    )

    # Get source
    source = Source()
    if source.is_unique(session, source_name):    
        source.name = source_name
        session.add(source)
        session.commit()
    sou = source.select_by_name(
        session,
        source_name
    )

    # Get experiment
    experiment = Experiment()
    if experiment.is_unique(
        session,
        experiment_type
    ):
        experiment.name = experiment_type
        session.add(experiment)
        session.commit()

    # Get samples
    for line in GUDglobals.parse_tsv_file(
        samples_file
    ):
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
        # Get sample
        samples.setdefault(cell_type, {
            "cell_or_tissue": cell_or_tissue,
            "add": add,
            "treatment": treatment,
            "cell_line": cell_line,
            "cancer": cancer
        })

    # For each sample...
    for s in samples:
        # Get sample
        sample = Sample()
        if sample.is_unique(
            session,
            samples[s]["cell_or_tissue"],
            samples[s]["treatment"],
            samples[s]["cell_line"],
            samples[s]["cancer"]
        ):
            sample.name =\
                samples[s]["cell_or_tissue"]
            sample.treatment =\
                samples[s]["treatment"]
            sample.cell_line =\
                samples[s]["cell_line"]
            sample.cancer =\
                samples[s]["cancer"]
            session.add(sample)
            session.commit()
        sam = sample.select_unique(
            session,
            samples[s]["cell_or_tissue"],
            samples[s]["treatment"],
            samples[s]["cell_line"],
            samples[s]["cancer"]
        )
        accession2sample.setdefault(
            s,
            sam.uid
        )
    print(accession2sample)
    exit(0)


    # Create table
    if not engine.has_table(table.__tablename__):
        # Create table
        table.__table__.create(bind=engine)

    # For each file...
    for bed_file in os.listdir(
        data_dir
    ):
        # Skip non-BED files
        if not bed_file.endswith(".bed") \
        and not bed_file.endswith(".bed.gz"):
            continue
        # Get TF name
        m = re.search(
            "^remap\d{4}_(.+)_all_.+.bed", bed_file)
        if m:
            # Initialize
            exp_dummy_dir = os.path.join(
                dummy_dir,
                "ChIP-seq.%s" % m.group(1)
            )
            # Skip if already done
            if os.path.isdir(exp_dummy_dir):
                continue
            # Create dummy dir
            if not os.path.isdir(exp_dummy_dir):
                os.mkdir(exp_dummy_dir)
            # For each line...
            for line in GUDglobals.parse_tsv_file(
                os.path.join(
                    data_dir, bed_file
                )
            ):
                # Initialize
                d, tf, s = line[3].split(".")
                tf = tf.split("_")
                s = s.split("_")
                # Fix sample
                if s[0] == "LNCAPABL":
                    s[0] == "LNCAP"
                # If sample exists...
                if s[0] in samples and tf[0] == m.group(1):
                    # Initialize
                    dummy_file = os.path.join(
                        exp_dummy_dir,
                        "%s.bed" % s[0]
                    )
                    # Write
                    GUDglobals.write(
                        dummy_file,
                        "\t".join(line)
                    )

    # For each directory...
    for directory in os.listdir(
        dummy_dir
    ):
        # Initialize
        dummy_file = os.path.join(
            dummy_dir,
            "%s.bed" % directory
        )
        exp_dummy_dir = os.path.join(
            dummy_dir, directory
        )
        # Skip if not directory
        if not os.path.isdir(
            exp_dummy_dir
        ): continue

        # Skip if completed
        completed_file = os.path.join(
            dummy_dir,
            "%s.completed" % directory
        )
        if not os.path.exists(
            completed_file
        ):
            # For each BED file...
            for bed_file in os.listdir(
                exp_dummy_dir
            ):
                # Initialize
                lines = []
                file_name = os.path.join(
                    exp_dummy_dir, bed_file
                )
                # Load BED
                a = BedTool(file_name)
                # For line in BED...
                for l in a:
                    # If BED12...
                    if len(l) >= 7:
                        lines.append(
                            "\t".join(
                                l[:7] + \
                                [
                                    str(int(l[6])+1),
                                    "0,0,0"
                                ]
                            )
                        )
                # Load BED from string
                a = BedTool(
                    "\n".join(lines), from_string=True
                )
                # Sort BED
                a = a.sort()
                # Save BED
                a.saveas(dummy_file)
                # Copy dummy BED to original
                shutil.copy(
                    dummy_file,
                    file_name
                )
                # Remove dummy BED
                os.remove(dummy_file)
            # Clean PyBedTools files
            cleanup(remove_all=True)
            # Completed
            open(completed_file, "a").close()

        # Cluster regions
        if cluster:
            # Initialize
            accession2sample = {}
            label2accession = {}
            regions = []
            bed_files = os.path.join(
                exp_dummy_dir,
                "files.txt"
            )
            table_file = os.path.join(
                exp_dummy_dir,
                "tableOfTables.txt"
            )
            cluster_file = os.path.join(
                exp_dummy_dir,
                "regCluster"
            )
            # Try clustering...
            try:
                # Create BED file list
                if not os.path.exists(bed_files):
                    # For each BED file...
                    for bed_file in os.listdir(
                        exp_dummy_dir
                    ):
                        # Skip non-BED files
                        if not bed_file.endswith(".bed"):
                            continue
                        # Skip regCluster BED file
                        if bed_file == "regCluster.bed":
                            continue
                        # Add file to list
                        GUDglobals.write(
                            bed_files,
                            os.path.join(
                                exp_dummy_dir,
                                bed_file
                            )
                        )
                # Make table of tables
                if not os.path.exists(table_file):
                    process = subprocess.check_output(
                        [
                            "regClusterMakeTableOfTables",
                            "uw01",
                            bed_files,
                            table_file
                        ],
                        stderr=subprocess.STDOUT
                    )
                # Make clusters
                if not os.path.exists(
                    "%s.cluster" % cluster_file
                ):
                    process = subprocess.check_output(
                        [
                            "regCluster",
                            table_file,
                            "%s.cluster" % cluster_file,
                            "%s.bed" % cluster_file
                        ],
                        stderr=subprocess.STDOUT
                    )
            except: pass
            exit(0)
            # For each accession, biosample...
            for accession, biosample in metadata[k]:
                # Get sample
                sample = Sample()
                if sample.is_unique(
                    session,
                    samples[biosample]["cell_or_tissue"],
                    samples[biosample]["treatment"],
                    samples[biosample]["cell_line"],
                    samples[biosample]["cancer"]
                ):
                    sample.name =\
                        samples[biosample]["cell_or_tissue"]
                    sample.treatment =\
                        samples[biosample]["treatment"]
                    sample.cell_line =\
                        samples[biosample]["cell_line"]
                    sample.cancer =\
                        samples[biosample]["cancer"]
                    session.add(sample)
                    session.commit()
                sam = sample.select_unique(
                    session,
                    samples[biosample]["cell_or_tissue"],
                    samples[biosample]["treatment"],
                    samples[biosample]["cell_line"],
                    samples[biosample]["cancer"]
                )
                accession2sample.setdefault(
                    accession,
                    sam.uid
                )
            # For each line...
            for line in GUDglobals.parse_tsv_file(
                table_file
            ):
                m = re.search("%s.%s/(\w+).bed" %\
                    (
                        experiment_type.replace(" ", "_"),
                        str(experiment_target)
                    ),
                    line[0]
                )
                if m:
                    label2accession.setdefault(
                        line[-1],
                        m.group(1)
                    )
            # For each line...
            for line in GUDglobals.parse_tsv_file(
                "%s.bed" % cluster_file
            ):
                # Get coordinates
                chrom = line[0]
                start = int(line[1])
                end = int(line[2])
                # Get region
                region = Region()
                if region.is_unique(
                    session,
                    chrom,
                    start,
                    end
                ):
                    # Insert region
                    region.bin = assign_bin(start, end)
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
            # For each line...
            for line in GUDglobals.parse_tsv_file(
                "%s.cluster" % cluster_file
            ):
                # Get region
                reg_uid = regions[int(line[0]) - 1] 
                # Get sample
                sam_uid = accession2sample[label2accession[line[-1]]]
                # Get accessibility feature
                if feat_type == "accessibility":
                    feat = DNAAccessibility()
                    is_unique = feat.is_unique(
                        session,
                        reg_uid,
                        sam_uid,
                        exp.uid,
                        sou.uid
                    )
                # Get histone feature
                if feat_type == "histone":
                    feat = HistoneModification()
                    is_unique = feat.is_unique(
                        session,
                        reg_uid,
                        sam_uid,
                        exp.uid,
                        sou.uid,
                        experiment_target
                    )
                # Get TF feature
                if feat_type == "tf":
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
                    feat.experimentID = exp.uid
                    feat.sourceID = sou.uid
                    if feat_type == "histone":
                        feat.histone_type = experiment_target
                    if feat_type == "tf":
                        feat.tf = experiment_target
                    session.add(feat)
                    session.commit()
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

#-------------#
# Main        #
#-------------#

if __name__ == "__main__": main()