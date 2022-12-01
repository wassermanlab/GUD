#!/usr/bin/env python

import argparse
from binning import assign_bin
import getpass
import os
import re
import sys
from sqlalchemy import create_engine
from sqlalchemy.orm import (
    scoped_session,
    sessionmaker
)
from sqlalchemy_utils import database_exists

# Import from GUD module
from GUD import GUDglobals
from GUD.ORM.chrom import Chrom
from GUD.ORM.dna_accessibility import DNAAccessibility
from GUD.ORM.enhancer import Enhancer
from GUD.ORM.experiment import Experiment
from GUD.ORM.histone_modification import HistoneModification
from GUD.ORM.region import Region
from GUD.ORM.sample import Sample
from GUD.ORM.source import Source
from GUD.ORM.tad import TAD
from GUD.ORM.tf_binding import TFBinding
#from .initialize import initialize_gud_db

usage_msg = """
usage: bed2gud.py --file [FILE ...] --feature STR
                  --experiment STR --sample STR --source STR 
                  [-h] [--histone STR | --enzyme STR | --tf STR]
                  [--cancer] [--cell-line] [--treatment]
                  [-d STR] [-H STR] [-p STR] [-P STR] [-u STR]
"""

help_msg = """%s

inserts genomic features from BED file into GUD. types of genomic
features include "accessibility", "enhancer", "histone", "tad" or
"tf".

  --file [FILE ...]   BED file(s)
  --feature STR       type of genomic feature (e.g. "enhancer")
  --experiment STR    experiment (e.g. "GRO-seq")
  --sample STR        sample (e.g. "B cell")
  --source STR        source (e.g. "PMID:29449408")

optional arguments:
  -h, --help          show this help message and exit
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

    # Mandatory args
    parser.add_argument("--file", nargs="*")
    parser.add_argument("--feature")
    parser.add_argument("--experiment")
    parser.add_argument("--sample")
    parser.add_argument("--source")

    # Optional args
    optional_group = parser.add_argument_group(
        "optional arguments"
    )
    optional_group.add_argument(
        "-h", "--help",
        action="store_true"
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

    # Print help
    if args.help:
        print(help_msg)
        exit(0)

    # Check mandatory arguments
    if (
        not args.file or \
        not args.feature or \
        not args.experiment or \
        not args.sample or \
        not args.source
    ):
        print(": "\
            .join(
                [
                    "%s\nbed2gud.py" % usage_msg,
                    "error",
                    "arguments \"--file\" \"--feature\" \"--experiment\" \"--sample\" \"--source\" are required\n"
                ]
            )
        )
        exit(0)

    # Check for invalid feature
    if args.feature not in feats:
        print(": "\
            .join(
                [
                    "%s\nbed2gud.py" % usage_msg,
                    "error",
                    "argument \"--feature\"",
                    "invalid choice",
                    "\"%s\" (choose from" % args.feature,
                    "%s)\n" % " "\
                    .join(["\"%s\"" % i for i in feats])
                ]
            )
        )
        exit(0)

    # Check feature type histone
    if args.feature == "histone":
        # Histone not provided!
        if not args.histone:
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
    if args.feature == "tad":
        # Enzyme not provided!
        if not args.enzyme:
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
    if args.feature == "tf":
        # TF not provided!
        if not args.tf:
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

    # Check MySQL password
    if not args.pwd:
        args.pwd = ""

def main():

    # Parse arguments
    args = parse_args()

    # For each BED file...
    for bed_file in args.file:
        # Insert BED file
        bed_to_gud_db(
            args.user,
            args.pwd,
            args.host,
            args.port,
            args.db,
            bed_file,
            args.feature,
            args.experiment,
            args.sample,
            args.source,
            args.histone,
            args.enzyme,
            args.tf,
            args.cancer,
            args.cell_line,
            args.treatment
        )

def bed_to_gud_db(user, pwd, host, port, db,
    bed_file, feat_type, experiment_type,
    sample_name, source_name, histone_type=None,
    restriction_enzyme=None, tf_name=None,
    cancer=False, cell_line=False,
    treatment=False):

    # Initialize
    features = []
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
    if feat_type == "tad":
        table = TAD()
    if feat_type == "tf":
        table = TFBinding()
    if not engine.has_table(table.__tablename__):
        # Create table
        table.__table__.create(bind=engine)

    # Get valid chromosomes
    chroms = Chrom.chrom_sizes(session)

    # For each line...
    for line in GUDglobals.parse_tsv_file(bed_file):
        # Skip if not enough elements
        if len(line) < 3: continue
        # Ignore non-standard chroms,
        # scaffolds, etc.
        if line[0] not in chroms:
            continue
        # Skip if not start or end
        if not line[1].isdigit(): continue
        if not line[2].isdigit(): continue
        # Initialize
        chrom = line[0]
        start = int(line[1])
        end = int(line[2])
        # If valid genomic coordinates...
        if start < end:
            # Get region
            region = Region()
            if region.is_unique(
                session,
                chrom,
                start,
                end
            ):
                # Insert region
                region.bin = assign_bin(
                    int(start),
                    int(end)
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
            # Get accessibility feature
            if feat_type == "accessibility":
                feat = DNAAccessibility()
                is_unique = feat.is_unique(
                    session,
                    reg.uid,
                    sam.uid,
                    exp.uid,
                    sou.uid
                )
            # Get enhancer feature
            if feat_type == "enhancer":
                feat = Enhancer()
                is_unique = feat.is_unique(
                    session,
                    reg.uid,
                    sam.uid,
                    exp.uid,
                    sou.uid
                )
            # Get histone feature
            if feat_type == "histone":
                feat = HistoneModification()
                is_unique = feat.is_unique(
                    session,
                    reg.uid,
                    sam.uid,
                    exp.uid,
                    sou.uid,
                    histone_type
                )
            # Get TAD feature
            if feat_type == "tad":
                feat = TAD()
                is_unique = feat.is_unique(
                    session,
                    reg.uid,
                    sam.uid,
                    exp.uid,
                    sou.uid,
                    restriction_enzyme
                )
            # Get TF feature
            if feat_type == "tf":
                feat = TFBinding()
                is_unique = feat.is_unique(
                    session,
                    reg.uid,
                    sam.uid,
                    exp.uid,
                    sou.uid,
                    tf_name
                )
            # Insert feature to GUD
            if is_unique:
                feat.regionID = reg.uid
                feat.sampleID = sam.uid
                feat.experimentID = exp.uid
                feat.sourceID = sou.uid
                if feat_type == "histone":
                    feat.histone_type = histone_type
                if feat_type == "tad":
                    feat.restriction_enzyme = restriction_enzyme
                if feat_type == "tf":
                    feat.tf = tf_name
                features.append(feat)

    session.add_all(features)
    session.commit()

#-------------#
# Main        #
#-------------#

if __name__ == "__main__": main()