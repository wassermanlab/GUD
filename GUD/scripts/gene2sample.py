#!/usr/bin/env python

import os, sys
import argparse
from Bio.SeqFeature import FeatureLocation
import ConfigParser
from sqlalchemy import create_engine, Index
from sqlalchemy.orm import Session
import shutil
import warnings

# Append GUD/OnTarget modules to path
module_path = os.path.join(os.path.dirname(__file__), os.pardir)
sys.path.append(os.path.join(module_path, os.pardir))
sys.path.append(os.path.join(module_path, "GUD"))

# Import from GUD/OnTarget
from GUD.ORM.tss import TSS
from OnTarget import OTglobals

# Read configuration file
config = ConfigParser.ConfigParser()
config_file = os.path.join(module_path, "config.ini")
config.read(config_file)

#-------------#
# Functions   #
#-------------#

def parse_args():
    """
    This function parses arguments provided via command
    line and returns an {argparse} object.
    """

    parser = argparse.ArgumentParser(description="describe what the script does...")

    parser.add_argument("--dummy-dir", default="/tmp/", help="Dummy directory (default = /tmp/)")

    gene_group = parser.add_mutually_exclusive_group(required=True)
    gene_group.add_argument("--gene", nargs="*", help="Gene symbol(s) (e.g. \"TH\")")
    gene_group.add_argument("--gene-file", help="File containing a list of gene names")

    # MySQL args
    mysql_group = parser.add_argument_group("exp. arguments")
    parser.add_argument("-a", "--all", action="store_true", help="Require exp. in all input sample(s) (default = False)")
    parser.add_argument("-g", "--group", action="store_true", help="Group exp. by gene (default = False)")
    parser.add_argument("-m", "--min-exp", type=float, default=10.0, help="Min. exp. in specified sample(s) (in TPM; default = 10.0)")
    parser.add_argument("-p", "--perc-exp", type=float, default=25.0,
        help="Min. percentage of exp. in input sample(s) (default = 25.0; i.e. expression in input sample(s) must account for at least 25 percent of total)")

    # MySQL args
    mysql_group = parser.add_argument_group("mysql arguments")
    mysql_group.add_argument("-d", "--db", default=config.get("MySQL", "db"),
        help="Database name (default = \"%s\")" % config.get("MySQL", "db"))
    mysql_group.add_argument("-H", "--host", default=config.get("MySQL", "host"),
        help="Host name (default = \"%s\")" % config.get("MySQL", "host"))
    mysql_group.add_argument("-P", "--port", default=config.get("MySQL", "port"),
        help="Port number (default = \"%s\")" % config.get("MySQL", "port"))
    mysql_group.add_argument("-u", "--user", default=config.get("MySQL", "user"),
        help="User name (default = \"%s\")" % config.get("MySQL", "user"))

    # Output args
    out_group = parser.add_argument_group("output arguments")
    out_group.add_argument("-o", "--out-file", help="Output file (default = stdout)")

    return parser.parse_args()

def select_differentially_expressed_tss(session, sample=[],
    min_tpm=0.0, perc_tpm=0.0, exp_in_all_samples=False,
    group_by_gene=False):
    """
    Query objects differentially expressed in samples.
    """

    # Initialize
    de_tss = {}
    all_tss = {}
    tss_samples = {}
    tss_expression = {}

    feats = TSS.select_by_sample(session, sample, min_tpm)

    # For each TSS...
    for tss in feats:
        # Add TSS sample
        tss_samples.setdefault((tss.gene, tss.tss), set())
        tss_samples[(tss.gene, tss.tss)].add(tss.cell_or_tissue)

    # Require expression in all samples
    if exp_in_all_samples:
        # For each TSS...
        for tss in frozenset(tss_samples):
            # For each sample...
            for s in sample:
                # If TSS not expressed in sample...
                if s not in tss_samples[tss]:
                    tss_samples.pop(tss, None)
                    break

    feats = TSS.select_by_multiple_tss(session, tss_samples.keys())

    # For each TSS...
    for tss in feats:
        # Add TSS
        tss_expression.setdefault((tss.gene, tss.tss), [0.0, 0.0])
        if tss.cell_or_tissue in sample:
            tss_expression[(tss.gene, tss.tss)][0] += tss.avg_tpm
        tss_expression[(tss.gene, tss.tss)][1] += tss.avg_tpm
        all_tss.setdefault((tss.gene, tss.tss), [])
        all_tss[(tss.gene, tss.tss)].append(tss)

    # For each TSS...
    for i in tss_expression:
        # Skip if not differentially expressed
        ptpm = tss_expression[i][0] * 100 / tss_expression[i][1]
        if ptpm < perc_tpm: continue
        # For each TSS...
        for tss in all_tss[i]:
            # If group by gene...
            if group_by_gene:
                de_tss.setdefault((tss.gene), [ptpm, {}, set()])
                de_tss[tss.gene][1].setdefault(tss.cell_or_tissue, [0.0, 0.0])
                de_tss[tss.gene][1][tss.cell_or_tissue][0] += tss.avg_tpm
                de_tss[tss.gene][1][tss.cell_or_tissue][1] = de_tss[tss.gene][1][tss.cell_or_tissue][0] * 100 / tss_expression[i][1]
                de_tss[tss.gene][2].add(tss.tss)
            # ... Else...
            else:
                de_tss.setdefault((tss.gene, tss.tss), [ptpm, {}])
                de_tss[(tss.gene, tss.tss)][1].setdefault(tss.cell_or_tissue, [tss.avg_tpm,
                    tss.avg_tpm * 100 / tss_expression[i][1]])

    return de_tss

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    # Parse arguments
    args = parse_args()

    # Initialize
    dummy_file = os.path.join(os.path.abspath(args.dummy_dir),
        "%s.%s.bed" % (os.path.basename(__file__), os.getpid()))

    # Fetch samples
    if args.sample:
        samples = args.sample
    elif args.sample_file:
        sample_file = os.path.abspath(args.sample_file)
        samples = [s for s in OTglobals.parse_file(sample_file)]
    else: samples = []

    # Establish a MySQL session
    db_name = "mysql://{}:@{}:{}/{}".format(args.user,
        args.host, args.port, args.db)
    try:
        engine = create_engine(db_name, echo=False)
        session = Session(engine)
    except:
        raise ValueError("Cannot connect to MySQL: %s" % db_name)

    # Get TSSs
    tss = select_differentially_expressed_tss(session, samples,
    args.min_exp, args.perc_exp, args.all, args.group)

    print(tss)
    exit(0)

    # For each gene...
    for gene in genes:
        # Get that gene's region
        region = get_gene_region(session, gene, samples,
            limit_by=args.limit)
        # Write
        OTglobals.write(dummy_file, region)

    # If output file...
    if args.out_file:
        shutil.copy(dummy_file, os.path.abspath(args.out_file))
    # ... Else, print on stdout
    else:
        for line in OTglobals.parse_file(dummy_file):
            OTglobals.write(None, line)

    # Delete dummy file
    os.remove(dummy_file)

  
