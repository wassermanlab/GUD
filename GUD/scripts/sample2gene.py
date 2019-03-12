#!/usr/bin/env python

import os, sys
import argparse
from Bio.SeqFeature import FeatureLocation
#import ConfigParser
import getpass
from sqlalchemy import create_engine, Index
from sqlalchemy.orm import Session
import shutil
import warnings

# Import from GUD module
from GUD2 import GUDglobals
from GUD2.ORM.gene import Gene
from GUD2.ORM.sample import Sample
from GUD2.ORM.tss import TSS

## Read configuration file
#config = ConfigParser.ConfigParser()
#config_file = os.path.join(module_path, "config.ini")
#config.read(config_file)

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

    sample_group = parser.add_mutually_exclusive_group(required=True)
    sample_group.add_argument("--sample", default=[], nargs="*", help="Sample(s) (limits RR search to using features from cells or tissues of the given sample; e.g. \"brain\")")
    sample_group.add_argument("--sample-file", help="File containing a list of samples")

    # MySQL args
    mysql_group = parser.add_argument_group("expression arguments")
    parser.add_argument("-a", "--all", action="store_true", help="require expression in all input sample(s) (default = False)")
    parser.add_argument("-g", "--group", action="store_true", help="group exp. by gene (default = False)")
    parser.add_argument("--min-exp", type=float, default=10.0,
        help="min. exp. in specified sample(s) (in TPM; default = 10.0)")
    parser.add_argument("--perc-exp", type=float, default=25.0,
        help="min. percentage of exp. in input sample(s) (default = 25.0; i.e. expression in input sample(s) must account for at least 25 percent of total)")

#    # MySQL args
#    mysql_group = parser.add_argument_group("mysql arguments")
#    mysql_group.add_argument("-d", "--db", default=config.get("MySQL", "db"),
#        help="Database name (default = \"%s\")" % config.get("MySQL", "db"))
#    mysql_group.add_argument("-H", "--host", default=config.get("MySQL", "host"),
#        help="Host name (default = \"%s\")" % config.get("MySQL", "host"))
#    mysql_group.add_argument("-P", "--port", default=config.get("MySQL", "port"),
#        help="Port number (default = \"%s\")" % config.get("MySQL", "port"))
#    mysql_group.add_argument("-u", "--user", default=config.get("MySQL", "user"),
#        help="User name (default = \"%s\")" % config.get("MySQL", "user"))

    # MySQL args
    mysql_group = parser.add_argument_group("mysql arguments")
    mysql_group.add_argument("-d", "--db", default="hg19",
        help="database name (default = hg19)")
    mysql_group.add_argument("-H", "--host", default="localhost",
        help="host name (default = localhost)")
    mysql_group.add_argument("-P", "--port", default=5506, type=int,
        help="port number (default = 5506)")

    user = getpass.getuser()
    mysql_group.add_argument("-u", "--user", default=user,
        help="user name (default = current user)")

    # Output args
    out_group = parser.add_argument_group("output arguments")
    out_group.add_argument("-o", "--out-file", help="output file (default = stdout)")

#    # Sample args
#    sample_group = parser.add_argument_group("sample arguments")
#    sample_group.add_argument("-c", "--cancer", action="store_true",
#        help="use \"cancer\" samples (default = False)")
#    sample_group.add_argument("-l", "--cell-line", action="store_true",
#        help="use samples from \"cell lines\" (default = False)")
#    sample_group.add_argument("-t", "--treatment", action="store_true",
#        help="use \"treated\" samples (default = False)")

    return parser.parse_args()

def main():

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
        samples = [s for s in GUDglobals.parse_file(sample_file)]
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
#    tss = select_differentially_expressed_tss(session, samples,
#        args.min_exp, args.perc_exp, args.all, args.group,
#        args.cancer, args.cell_line, args.treatment)
    tss = select_differentially_expressed_tss(session, samples,
        args.min_exp, args.perc_exp, args.all, args.group)

    # For each TSS...
    for i in sorted(tss, key=lambda x: tss[x][0], reverse=True):
        print(i, tss[i][0])
#        # If group by gene...
#        if args.group:

#        # Write
#        OTglobals.write(dummy_file, region)
#
#    # If output file...
#    if args.out_file:
#        shutil.copy(dummy_file, os.path.abspath(args.out_file))
#    # ... Else, print on stdout
#    else:
#        for line in OTglobals.parse_file(dummy_file):
#            OTglobals.write(None, line)
#
#    # Delete dummy file
#    os.remove(dummy_file)

#def select_differentially_expressed_tss(session, names=[],
#    min_tpm=0.0, perc_tpm=0.0, exp_in_all_samples=False,
#    group_by_gene=False, cancer=False, cell_line=False,
#    treatment=False):
def select_differentially_expressed_tss(session, names=[],
    min_tpm=0.0, perc_tpm=0.0, exp_in_all_samples=False,
    group_by_gene=False):
    """Select TSSs differentially expressed in samples.
    """

    # Initialize
    tss_in_samples = {}
    diff_exp_tss = {}

    # Get genes
    all_genes = set(Gene.get_all_gene_symbols(session))

    # Get samples
    feats = Sample.select_by_names(session, names)
    samples = [f.uid for f in feats]
    feats = Sample.select_by_exp_conditions(session,
        cancer, cell_line, treatment)
    all_samples = [f.uid for f in feats]

    # Get TSSs
    feats = TSS.select_by_samples(session, samples, min_tpm)

    # For each TSS...
    for tss in feats:
        # Exclude non gene TSSs
        if tss.gene not in all_genes: continue
        # Add TSS sample
        tss_in_samples.setdefault((tss.gene, tss.tss), set())
        tss_in_samples[(tss.gene, tss.tss)].add(tss.sampleID)

    print(len(tss_in_samples))
    # Require expression in all samples
    if exp_in_all_samples:
        # Get TSS samples
        all_tss_samples = TSS.get_all_samples(session)
        # Get valid samples (i.e. available in TSS table)
        valid_samples = set(samples).intersection(
            set(all_tss_samples))
        # For each TSS...
        for tss in frozenset(tss_in_samples):
            # If TSS not expressed in sample...
            if len(tss_in_samples[tss]) < len(valid_samples):
                # Remove TSS
                tss_in_samples.pop(tss, None)

    # For each TSS...
    for gene, tss in sorted(tss_in_samples):
        # Get TSSs
        feats = TSS.select_by_tss(
            session, gene, tss, all_samples)
        # Get relative exp.
        tpm_in_samples = sum([f.avg_tpm for f in feats if f.sampleID in tss_in_samples[(gene, tss)]])
        total_tpm = sum([f.avg_tpm for f in feats])
        # If differentially exp.
        rel_tpm = tpm_in_samples * 100.0 / total_tpm
        if rel_tpm >= perc_tpm:
            diff_exp_tss.setdefault((gene, tss), [rel_tpm, feats])

    for tss in sorted(diff_exp_tss, key=lambda x: diff_exp_tss[x][0], reverse=True):
        print(tss, diff_exp_tss[tss][0])
    exit(0)
#    feats = TSS.select_by_multiple_tss(session,
#        tss_samples.keys(), all_samples)
    
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

    main()