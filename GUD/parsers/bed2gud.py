#!/usr/bin/env python

import argparse
from binning import assign_bin
import getpass
import pybedtools
import os
import re
import sys
from sqlalchemy import create_engine
from sqlalchemy.orm import (
    scoped_session,
    sessionmaker
)
from sqlalchemy_utils import database_exists
import warnings

# Import from GUD module
from GUD import GUDglobals
from GUD.ORM.dna_accessibility import DNAAccessibility
from GUD.ORM.enhancer import Enhancer
from GUD.ORM.experiment import Experiment
from GUD.ORM.histone_modification import HistoneModification
from GUD.ORM.region import Region
from GUD.ORM.sample import Sample
from GUD.ORM.source import Source
from GUD.ORM.tad import TAD
from GUD.ORM.tf_binding import TFBinding

usage_msg = """
usage: bed2gud.py file feat_type experiment sample source
                  [-h] [-c] [--dummy-dir DIR]
                  [[--histone STR] | [--enzyme STR] | [--tf STR]]
                  [-d STR] [-H STR] [-p STR] [-P STR] [-u STR]
"""

help_msg = """%s

inserts genomic features from BED file into GUD. types of genomic
features include "accessibility", "enhancer", "histone", "tad" or
"tf".

  file                BED file(s)
  feat_type           type of genomic feature
  experiment          genomic experiment (e.g. "GRO-seq")
  sample              sample (e.g. "brain")
  source              source name (e.g. "PMID:29449408")

optional arguments:
  -h, --help          show this help message and exit
  -c, --cluster       cluster genome regions by UCSC's regCluster
                      (default = False)
  --dummy-dir DIR     dummy directory (default = "/tmp/")
  --histone STR       histone type (e.g. "H3K27ac")
  --enzyme STR        restriction enzyme (e.g. "HindIII")
  --tf STR            TF name (e.g. "FOS")

sample arguments:
  --cancer            label sample as cancer (default = False)
  --cell-line         label sample as cell line (default = False)
  --treatment         label sample as treated (default = False)

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
    optional_group.add_argument("--histone")
    optional_group.add_argument("--enzyme")
    optional_group.add_argument("--tf")

    # Sample args
    sample_group = parser.add_argument_group(
        "sample arguments"
    )
    sample_group.add_argument(
        "--cancer",
        action="store_true"
    )
    sample_group.add_argument(
        "--cell-line",
        action="store_true"
    )
    sample_group.add_argument(
        "--treatment",
        action="store_true"
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

    if "-h" in sys.argv or "--help" in sys.argv:
        print(help_msg)
        exit(0)

    # Mandatory arguments
    parser.add_argument("file", nargs="*")
    parser.add_argument("feat_type")
    parser.add_argument("experiment")
    parser.add_argument("sample")
    parser.add_argument("source")

    args = parser.parse_args()

    check_args(args)

    return args

def check_args(args):
    """
    This function checks an {argparse} object.
    """

    # Initialize
    feats = [
        "accessibility",
        "enhancer",
        "histone",
        "tad",
        "tf"
    ]

    # Check for invalid feature
    if args.feat_type not in feats:
        print(": "\
            .join(
                [
                    "%s\nbed2gud.py" % usage_msg,
                    "error",
                    "argument \"feat_type\"",
                    "invalid choice",
                    "\"%s\" (choose from" % args.feat_type,
                    "%s)\n" % ", "\
                    .join(["\"%s\"" % i for i in feats])
                ]
            )
        )
        exit(0)

    # Check feature type histone
    if args.feat_type == "histone":
        if args.histone: continue
        # Histone not provided!
        print(": "\
            .join(
                [
                    "%s\nbed2gud.py" % usage_msg,
                    "error",
                    "argument \"--histone\"",
                    "missing argument\n"
                ]
            )
        )
        exit(0)

    # Check feature type TAD
    if args.feat_type == "tad":
        if args.enzyme: continue
        # Enzyme not provided!
        print(": "\
            .join(
                [
                    "%s\nbed2gud.py" % usage_msg,
                    "error",
                    "argument \"--enzyme\"",
                    "missing argument\n"
                ]
            )
        )
        exit(0)

    # Check feature type TF
    if args.feat_type == "tf":
        if args.tf: continue
        # TF not provided!
        print(": "\
            .join(
                [
                    "%s\nbed2gud.py" % usage_msg,
                    "error",
                    "argument \"--tf\"",
                    "missing argument\n"
                ]
            )
        )
        exit(0)

def main():

    # Parse arguments
    args = parse_args()

    # Insert BED file to GUD database
    insert_bed_to_gud_db(
        args.user,
        args.pwd,
        args.host,
        args.port,
        args.db,
        args.file,
        args.feat_type,
        args.experiment,
        args.sample,
        args.source,
        args.histone,
        args.enzyme,
        args.tf
    )

def insert_bed_to_gud_db(user, pwd, host, port,
    db, bed_files, feat_type, experiment_type,
    sample_name, source_name, cluster=False,
    dummy_dir="/tmp/" histone_type=None,
    restriction_enzyme=None, tf_name=None,
    cancer=False, cell_line=False,
    treatment=False):

    # Initialize
    lines = []
    db_name = "mysql://{}:{}@{}:{}/{}".format(
        user, pwd, host, port, db
    )
    if not database_exists(db_name):
        raise ValueError("GUD database does not exist!!!\n\t%s" % db_name)
    session = scoped_session(sessionmaker())
    engine = create_engine(
        db_name,
        echo=False
    )
    session.remove()
    session.configure(
        bind=engine,
        autoflush=False,
        expire_on_commit=False
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
    exp = experiment.select_by_name(
        session,
        experiment_type
    )

    # Get sample
    sample = Sample()
    if sample.is_unique(
        session,
        sample_name,
        int(treatment),
        int(cell_line),
        int(cancer)
    ):
        sample.name = sample_name
        sample.treatment = int(treatment)
        sample.cell_line = int(cell_line)
        sample.cancer = int(cancer)
        session.add(sample)
        session.commit()
    sam = sample.select_unique(
        session,
        sample_name,
        int(treatment),
        int(cell_line),
        int(cancer)
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

    # Create table
    if feat_type == "accessibility":
        table = DNAAccessibility()
    if feat_type == "enhancer":
        table = Enhancer()
    if feat_type == "histone":
        table = HistoneModification()
    if feat_type == "TAD":
        table = TFBinding()
    if feat_type == "tf":
        table = TFBinding()
    if not engine.has_table(table.__tablename__):
        # Create table
        table.__table__.create(bind=engine)

    # For each BED file...
    for bed_file in bed_files:
        # For each line...
        for line in GUDglobals.parse_tsv_file(bed_file):
            # Skip if not enough elements
            if len(line) < 3: continue
            # Ignore non-standard chroms, scaffolds, etc.
            m = re.search("^chr(\S+)$", line[0])
            if not m.group(1) in GUDglobals.chroms:
                continue
            # Skip if not start or end
            if not line[1].isdigit(): continue
            if not line[2].isdigit(): continue
            # If valid genomic coordinates...
            if int(line[1]) < int(line[2]):
                lines.append("\t".join(line[:3]))

    # If lines...
    if lines:
        # Create BED object
        bed_obj = pybedtools.BedTool(
            "\n".join(lines),
            from_string=True
        )
        
        # Cluster...
        if cluster:
            # 2 be implemented
            # copy from encode2gud
            pass

        # Do not cluster...
        else:
            # Sort BED object
            for chrom, start, end in bed_obj.sort():
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

                if feat_type == "accessibility":
                    pass

                if feat_type == "enhancer":
                    pass

                if feat_type == "histone":
                    pass

                if feat_type == "tad":
                    pass

                if feat_type == "tf":
                    pass

#-------------#
# Main        #
#-------------#

if __name__ == "__main__": main()