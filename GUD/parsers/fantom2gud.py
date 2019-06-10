#!/usr/bin/env python

import argparse
from binning import assign_bin
import getpass
import re
from sqlalchemy import create_engine
from sqlalchemy.orm import (
    scoped_session,
    sessionmaker
)
from sqlalchemy_utils import database_exists
try: from urllib2 import unquote
except: from urllib.parse import unquote

# Import from GUD module
from GUD import GUDglobals
from GUD.ORM.enhancer import Enhancer
from GUD.ORM.experiment import Experiment
from GUD.ORM.expression import Expression
from GUD.ORM.region import Region
from GUD.ORM.sample import Sample
from GUD.ORM.source import Source
from GUD.ORM.tss import TSS

usage_msg = """
usage: %s  --matrix FILE --samples FILE --feature STR
                      [-h] [--source STR]
                      [-d STR] [-H STR] [-p STR] [-P STR] [-u STR]
""" % \
os.path.basename(__file__)

help_msg = """%s
inserts "enhancer" and "tss" features from the FANTOM5 consortium
into GUD. "--matrix" refers to:

    1) "hg19_permissive_enhancers_expression_rle_tpm.csv.gz"; and
    2) "hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt.gz"

  --matrix FILE       expression matrix TPM/RLE normalized across
                      all FANTOM libraries
  --samples FILE      FANTOM samples (manually-curated)
  --feature STR       type of genomic feature

optional arguments:
  -h, --help          show this help message and exit
  --source STR        source name (default = "FANTOM5")

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
    parser.add_argument("--matrix")
    parser.add_argument("--samples")
    parser.add_argument("--feature")

    # Optional args
    optional_group = parser.add_argument_group(
        "optional arguments"
    )
    optional_group.add_argument(
        "-h", "--help",
        action="store_true"
    )
    optional_group.add_argument(
        "--source",
        default="FANTOM5"
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
        "enhancer",
        "tss"
    ]

    # Print help
    if args.help:
        print(help_msg)
        exit(0)

    # Check mandatory arguments
    if (
        not args.matrix or \
        not args.samples or \
        not args.feature
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
                    "arguments \"--matrix\" \"--samples\" \"--feature\" are required\n"
                ]
            )
        )
        exit(0)

    # Check for invalid feature
    if args.feature not in feats:
        print(": "\
            .join(
                [
                    "%s\n%s" % \
                        (
                            usage_msg,
                            os.path.basename(__file__)
                        ),
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

    # Check MySQL password
    if not args.pwd:
        args.pwd = ""

def main():

    # Parse arguments
    args = parse_args()

    # Insert FANTOM data
    fantom_to_gud_db(
        args.user,
        args.pwd,
        args.host,
        args.port,
        args.db,
        args.matrix,
        args.samples,
        args.feature,
        args.source
    )

def fantom_to_gud_db(user, pwd, host, port, db,
    matrix_file, samples_file, feat_type,
    source_name="FANTOM5"):

    # Initialize
    samples = {}
    sample_ids = []
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

    # Get samples
    for line in GUDglobals.parse_tsv_file(
        samples_file
    ):
        # Initialize 
        sample_id = line[0]
        sample_name = line[1]
        cell_or_tissue = line[2]
        if line[3] == "Yes": add = True
        else: add = False
        if line[4] == "Yes": treatment = True
        else: treatment = False
        if line[5] == "Yes": cell_line = True
        else: cell_line = False
        if line[6] == "Yes": cancer = True
        else: cancer = False
        # Get sample
        samples.setdefault(sample_id,
            {
                "sample_name": sample_name,
                "cell_or_tissue": cell_or_tissue,
                "add": add,
                "treatment": treatment,
                "cell_line": cell_line,
                "cancer": cancer
            }
        )

    # Get experiment
    experiment = Experiment()
    if experiment.is_unique(session, "CAGE"):    
        experiment.name = "CAGE"
        session.add(experiment)
        session.commit()
    exp = experiment.select_by_name(
        session,
        "CAGE"
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
    if feat_type == "enhancer":
        tpms_start_at = 1
        table = Enhancer()
        if not engine.has_table(
            table.__tablename__
        ):
            # Create table
            table.__table__.create(
                bind=engine
            )
        lines = GUDglobals.parse_csv_file(
            matrix_file
        )

    if feat_type == "tss":
        tpms_start_at = 7
        table = TSS()
        if not engine.has_table(
            table.__tablename__
        ):
            # Create table
            table.__table__.create(
                bind=engine
            )
        table = Expression()
        if not engine.has_table(
            table.__tablename__
        ):
            # Create table
            table.__table__.create(
                bind=engine
            )
        lines = GUDglobals.parse_tsv_file(
            matrix_file
        )

    # For each line...
    for line in lines:
        # Skip comments
        if line[0].startswith("#"): continue
        # Get sample IDs
        if len(sample_ids) == 0:
            for sample in line[tpms_start_at:]:
                m = re.search(
                    "(CNhs\d+)",
                    unquote(sample)
                )
                sample_ids.append(m.group(1))
        # If enhancer/TSS...
        elif line[0].startswith("chr") or\
             line[0].startswith("\"chr"):
            # Initialize
            data = {}
            features = []
            sampleIDs = []
            avg_expression_levels = []
            # Get coordinates
            if feat_type == "enhancer":
                m = re.search(
                    "(chr\S+)\:(\d+)\-(\d+)",
                    line[0]
                )
                chrom = m.group(1)
                start = int(m.group(2))
                end = int(m.group(3))
                strand = None
            if feat_type == "tss":
                # Initialize
                peak_ids = set()
                m = re.search(
                    "(chr\S+)\:(\d+)\.\.(\d+),(\S)",
                    line[0]
                )
                chrom = m.group(1)
                start = int(m.group(2))
                end = int(m.group(3))
                strand = m.group(4)
                for peak in line[1].split(","):
                    m = re.search(
                        "p(\d+)@(\w+)",
                        peak
                    )
                    if m:
                        peak_ids.add((
                            m.group(2), m.group(1)
                        ))
                    else:
                        peak_ids.add((None, 1))
            # Ignore non-standard chroms,
            # scaffolds, etc.
            m = re.search("^chr(\S+)$", chrom)
            if not m.group(1) in GUDglobals.chroms:
                continue
            # Get region
            region = Region()
            if region.is_unique(
                session,
                chrom,
                start,
                end,
                strand
            ):
                # Insert region
                region.bin = assign_bin(start, end)
                region.chrom = chrom
                region.start = start
                region.end = end
                region.strand = strand
                session.add(region)
                session.commit()
            reg = region.select_unique(
                session,
                chrom,
                start,
                end,
                strand
            )
            # For each sample...
            for i in range(tpms_start_at, len(line)):
                # Skip sample
                sample_id = sample_ids[
                    i - tpms_start_at
                ]
                if not samples[sample_id]["add"]:
                    continue
                # Get data
                name =\
                    samples[sample_id]["cell_or_tissue"]
                treatment =\
                    samples[sample_id]["treatment"]
                cell_line =\
                    samples[sample_id]["cell_line"]
                cancer =\
                    samples[sample_id]["cancer"]
                data.setdefault((
                        name,
                        treatment,
                        cell_line,
                        cancer
                    ), []
                )
                data[(
                        name,
                        treatment,
                        cell_line,
                        cancer
                    )
                ].append(float(line[i]))
            # For each sample...
            for s in data:
                # Initialize
                name = s[0]
                treatment = s[1]
                cell_line = s[2]
                cancer = s[3]
                # Get sample
                sample = Sample()
                if sample.is_unique(
                    session,
                    name,
                    treatment,
                    cell_line,
                    cancer
                ):
                    sample.name = name
                    sample.treatment = treatment
                    sample.cell_line = cell_line
                    sample.cancer = cancer
                    session.add(sample)
                    session.commit()
                sam = sample.select_unique(
                    session,
                    name,
                    treatment,
                    cell_line,
                    cancer
                )
                # Skip if feature not expressed in sample
                avg_expression_level = float(
                    sum(data[
                        name,
                        treatment,
                        cell_line,
                        cancer
                    ]) / len(data[
                        name,
                        treatment,
                        cell_line,
                        cancer
                    ])
                )
                if avg_expression_level > 0:
                    sampleIDs.append(sam.uid)
                    avg_expression_levels.append(
                        "%.3f" % avg_expression_level
                    )
            # Get TSS
            if feat_type == "tss":
                # Initialize
                tss_ids = set()
                for gene, tss_id in peak_ids:
                    tss = TSS()
                    if tss.is_unique(
                        session,
                        reg.uid,
                        sou.uid,
                        exp.uid,
                        gene,
                        tss_id
                    ):
                        tss.regionID = reg.uid
                        tss.sourceID = sou.uid
                        tss.experimentID = exp.uid
                        tss.gene = gene
                        tss.tss = tss_id
                        tss.sampleIDs = "{},".format(
                            ",".join(map(str, sampleIDs))
                        )
                        tss.avg_expression_levels = "{},".format(
                            ",".join(avg_expression_levels)
                        )
                        session.add(tss)
                        session.commit()
                    tss = tss.select_unique(
                        session,
                        reg.uid,
                        sou.uid,
                        exp.uid,
                        gene,
                        tss_id
                    )
                    tss_ids.add(tss)
            # For each sample...
            for i in range(len(sampleIDs)):
                if feat_type == "enhancer":
                    enhancer = Enhancer()
                    if enhancer.is_unique(
                        session,
                        reg.uid,
                        sou.uid,
                        sampleIDs[i],
                        exp.uid
                    ):
                        enhancer.regionID = reg.uid
                        enhancer.sourceID = sou.uid
                        enhancer.sampleID = sampleIDs[i]
                        enhancer.experimentID = exp.uid
                        # Append to features
                        features.append(enhancer)
                if feat_type == "tss":
                    for tss in tss_ids:
                        expression = Expression()
                        if expression.is_unique(
                            session,
                            tss.uid,
                            sampleIDs[i]
                        ):
                            expression.tssID = tss.uid
                            expression.sampleID = sampleIDs[i]
                            expression.avg_expression_level =\
                                avg_expression_levels[i]
                            # Append to features
                            features.append(expression)
            # Insert features
            session.add_all(features)
            session.commit()

#-------------#
# Main        #
#-------------#

if __name__ == "__main__": main()