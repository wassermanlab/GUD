#!/usr/bin/env python

import argparse
import os
import shutil

# Import from GUD
from GUD import GUDglobals
from GUD.ORM.chrom import Chrom
from GUD.ORM.expression import Expression
from GUD.ORM.gene import Gene
from GUD.ORM.seq_feature import SeqFeature
from GUD.ORM.tss import TSS

#-------------#
# Functions   #
#-------------#

def parse_args():
    """
    This function parses arguments provided via command
    line and returns an {argparse} object.
    """

    parser = argparse.ArgumentParser(description="identifies gene(s) differentially expressed in sample(s).")

    parser.add_argument("--dummy-dir", default="/tmp/",
        help="dummy directory (default = /tmp/)")

    sample_group = parser.add_mutually_exclusive_group(required=True)
    sample_group.add_argument("--sample", default=[], nargs="*",
        help="sample(s) (e.g. \"brain\")")
    sample_group.add_argument("--sample-file",
        help="file containing a list of samples")

    # MySQL args
    mysql_group = parser.add_argument_group("exp. arguments")
    parser.add_argument("-a", "--all", action="store_true",
        help="require expression in all input sample(s) (default = False)")
    parser.add_argument("-g", "--group", action="store_true",
        help="group expression by gene (default = False)")
    parser.add_argument("-m", "--min-exp", type=float, default=GUDglobals.min_exp,
        help="min. expression in input sample(s) (in TPM; default = %s)" % GUDglobals.min_exp)
    parser.add_argument("-p", "--min-percent-exp", type=float, default=GUDglobals.min_percent_exp,
        help="min. percentage of expression in input sample(s) (default = %s)" % GUDglobals.min_percent_exp)

    # MySQL args
    mysql_group = parser.add_argument_group("mysql arguments")
    mysql_group.add_argument("-d", "--db", default=GUDglobals.db_name,
        help="database name (default = \"%s\")" % GUDglobals.db_name)
    mysql_group.add_argument("-H", "--host", default=GUDglobals.db_host,
        help="host name (default = \"%s\")" % GUDglobals.db_host)
    mysql_group.add_argument("-p", "--passwd", metavar="PASS",
        help="password (default = ignore this option)")
    mysql_group.add_argument("-P", "--port", default=GUDglobals.db_port,
        help="port number (default = \"%s\")" % GUDglobals.db_port)
    mysql_group.add_argument("-u", "--user", default=GUDglobals.db_user,
        help="user name (default = \"%s\")" % GUDglobals.db_user)

    # Output args
    out_group = parser.add_argument_group("output arguments")
    out_group.add_argument("-o", "--out-file",
        help="output file (default = stdout)")

    args = parser.parse_args()

    return args

def main():

    # Parse arguments
    args = parse_args()

    # Initialize
    dummy_file = os.path.join(
        os.path.abspath(args.dummy_dir),
        "%s.%s.txt" % (os.path.basename(__file__),
        os.getpid()))

    # Fetch samples
    samples = []
    if args.sample:
        samples = args.sample
    elif args.sample_file:
        sample_file = os.path.abspath(args.sample_file)
        for s in GUDglobals.parse_file(sample_file):
            samples.append(s)
    else:
        raise ValueError("No sample(s) was provided!")

    print(samples)
    exit(0)

    # Establish SQLalchemy session with GUD
    session = GUDglobals.establish_GUD_session(
        args.user,
        args.passwd,
        args.host,
        args.port,
        args.db
    )

    # Get TSSs
    tss = get_differentially_expressed_tss(
        session,
        samples,
        args.min_exp,
        args.perc_exp,
        args.all,
        args.group
    )

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

def get_differentially_expressed_tss(session,
    sample=[], min_tpm=0.0, perc_tpm=0.0,
    exp_in_all_samples=False, group_by_gene=False):
    """
    Identifies TSSs differentially expressed in samples.
    """

    # Initialize
    de_tss = {}
    all_tss = {}
    tss_samples = {}
    tss_expression = {}

    feats = TSS.select_by_sample(session, sample, min_tpm)
    print(len(feats))

    gene_symbols = set(Gene.select_all_gene_symbols(session))

    # For each TSS...
    for tss in feats:
        # Exclude non gene TSSs
        if tss.gene not in gene_symbols: continue
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
    print(len(feats))

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

if __name__ == "__main__": main()