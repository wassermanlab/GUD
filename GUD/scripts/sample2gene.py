#!/usr/bin/env python

import argparse
import os
import re
import shutil
import sys

# Import from GUD
from GUD import GUDglobals
from GUD.ORM.expression import Expression
from GUD.ORM.gene import Gene
from GUD.ORM.sample import Sample
from GUD.ORM.tss import TSS

usage_msg = """
usage: sample2gene.py (--sample [STR ...] | --sample-file FILE)
                      [-h] [--dummy-dir DIR] [-g] [-o FILE]
                      [-a] [--percent FLT] [--tpm FLT] [--tss INT]
                      [-d STR] [-H STR] [-p STR] [-P STR] [-u STR]
"""

help_msg = """%s
identifies one or more genes selectively expressed in samples.

  --sample [STR ...]  sample(s) (e.g. "B cell")
  --sample-file FILE  file containing a list of samples

optional arguments:
  -h, --help          show this help message and exit
  --dummy-dir DIR     dummy directory (default = "/tmp/")
  -g, --group         group by gene (if True, return the most
                      selectively expressed TSS of each gene;
                      default = False)
  -o FILE             output file (default = stdout)

expression arguments:
  -a, --all           expression in all samples (default = False)
  --percent FLT       min. percentile of expression for TSS in
                      input samples (default = %s)
  --tpm FLT           min. expression levels (in TPM) for TSS in
                      input samples (default = %s)
  --tss INT           max. number of TSSs to return (if 0, return
                      all TSSs; default = %s)

mysql arguments:
  -d STR, --db STR    database name (default = "%s")
  -H STR, --host STR  host name (default = "%s")
  -p STR, --pwd STR   password (default = ignore this option)
  -P STR, --port STR  port number (default = %s)
  -u STR, --user STR  user name (default = "%s")
""" % \
(
    usage_msg,
    GUDglobals.min_percent,
    GUDglobals.min_tpm,
    GUDglobals.max_tss,
    GUDglobals.db_name,
    GUDglobals.db_host,
    GUDglobals.db_port,
    GUDglobals.db_user
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

    # Mandatory arguments
    parser.add_argument(
        "--sample", nargs="*"
    )
    parser.add_argument("--sample-file")

    # Optional args
    optional_group = parser.add_argument_group(
        "optional arguments"
    )
    optional_group.add_argument(
        "-h", "--help",
        action="store_true"
    )
    optional_group.add_argument(
        "--dummy-dir",
        default="/tmp/"
    )
    optional_group.add_argument(
        "-g", "--group",
        action="store_true"
    )
    optional_group.add_argument("-o")

    # Expression args
    exp_group = parser.add_argument_group(
        "expression arguments"
    )
    exp_group.add_argument(
        "-a", "--all",
        action="store_true"
    )
    exp_group.add_argument(
        "--percent",
        default=GUDglobals.min_percent
    )
    exp_group.add_argument(
        "--tpm",
        default=GUDglobals.min_tpm
    )
    exp_group.add_argument(
        "--tss",
        default=GUDglobals.max_tss
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
        default=GUDglobals.db_host
    )
    mysql_group.add_argument("-p", "--pwd")
    mysql_group.add_argument(
        "-P", "--port",
        default=GUDglobals.db_port
    )
    mysql_group.add_argument(
        "-u", "--user",
        default=GUDglobals.db_user
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
        not args.sample and \
        not args.sample_file
    ):
        print(": "\
            .join(
                [
                    "%s\nsample2gene.py" % usage_msg,
                    "error",
                    "one of the arguments \"--sample\" \"--sample-file\" is required\n"
                ]
            )
        )
        exit(0)

    if args.sample and args.sample_file:
        print(": "\
            .join(
                [
                    "%s\nsample2gene.py" % usage_msg,
                    "error",
                    "arguments \"--sample\" \"--sample-file\"",
                    "expected one argument\n"
                ]
            )
        )
        exit(0)

    # Check for invalid percent
    try:
        args.percent = float(args.percent)
    except:
        print(": "\
            .join(
                [
                    "%s\nsample2gene.py" % usage_msg,
                    "error",
                    "argument \"--percent\"",
                    "invalid float value",
                    "\"%s\"\n" % args.percent
                ]
            )
        )
        exit(0)
        
    # Check for invalid TPM
    try:
        args.tpm = float(args.tpm)
    except:
        print(": "\
            .join(
                [
                    "%s\nsample2gene.py" % usage_msg,
                    "error",
                    "argument \"--tpm\"",
                    "invalid float value",
                    "\"%s\"\n" % args.tpm
                ]
            )
        )
        exit(0)

    # Check for invalid TSS
    try:
        args.tss = int(args.tss)
    except:
        print(": "\
            .join(
                [
                    "%s\nsample2gene.py" % usage_msg,
                    "error",
                    "argument \"--tss\"",
                    "invalid int value",
                    "\"%s\"\n" % args.tss
                ]
            )
        )
        exit(0)

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
        for sample in GUDglobals.parse_file(sample_file):
            samples.append(sample)
    else:
        raise ValueError(
            "No samples were provided!!!"
        )

    # Establish SQLalchemy session with GUD
    session = GUDglobals.establish_GUD_session(
        args.user,
        args.pwd,
        args.host,
        args.port,
        args.db
    )

    # Get selectively expressed TSSs
    tss = get_selectively_expressed_tss(
        session,
        samples,
        args.tpm,
        args.percent,
        args.all
    )

    # If selectively expressed TSSs...
    if tss:

        # Initialize
        genes = set()
        tss_count = 0
        if args.tss <= 0:
            args.tss = len(tss)

        # For each TSS...
        for t in tss:
            # Skip if enough TSSs
            if tss_count == args.tss:
                break
            # If group by gene...
            if args.group:
                # Skip if there is a better
                # TSS for this gene
                if t.qualifiers["gene"] in \
                    genes: continue
            # Write
            GUDglobals.write(dummy_file, t)
            # Add gene and increase count
            genes.add(t.qualifiers["gene"])
            tss_count += 1

        # If output file...
        if args.o:
            shutil.copy(
                dummy_file,
                os.path.abspath(args.o)
            )
        # ... Else, print on stdout...
        else:
            for line in GUDglobals.parse_file(
                dummy_file
            ):
                GUDglobals.write(None, line)

        # Delete dummy file
        os.remove(dummy_file)

    else:
        raise ValueError(
            "No selectively expressed genes found!!!"
        )

def get_selectively_expressed_tss(session,
    samples=[], tpm_exp=100, percent_exp=0.0,
    exp_in_all_samples=False):
    """
    Identifies TSSs from genes selectively
    expressed in samples.
    """

    # Initialize
    name2uid = {}
    uid2name = {}
    tssIDs = set()
    sample2tss = []
    sel_exp_tss = []
    isfloat = re.compile("\d+(\.\d+)?")

    # For each sample...
    for feat in Sample.select_by_names(session):
        uid2name.setdefault(feat.uid, feat.name)
        name2uid.setdefault(feat.name, feat.uid)

    # Get genes
    genes = set(Gene.get_all_gene_symbols(session))

    # For each sample...
    for sample in samples:
        sample2tss.append(set())
    # For each TSS...
    for tss in Expression.select_by_samples(
        session, samples, tpm_exp
    ):
        # Get index
        i = samples.index(tss.Sample.name)
        sample2tss[i].add(tss.Expression.tssID)
        tssIDs.add(tss.Expression.tssID)
        # Enables search for cell line
        # specific TSSs
        uid2name.setdefault(
            tss.Sample.uid, tss.Sample.name)
        name2uid.setdefault(
            tss.Sample.name, tss.Sample.uid)
    # If required expression in all samples...
    if exp_in_all_samples:
        tssIDs = list(
            set.intersection(*sample2tss)
        )
    # ... Else...
    else:
        tssIDs = list(tssIDs)

    # If there are TSS IDs...
    if tssIDs:
        # For each TSS...
        for tss in TSS.select_by_uids(session,
            tssIDs, as_genomic_feature=True):
            # If gene TSS...
            if tss.qualifiers["gene"] in genes:
                # Initialize
                expression = {}
                bg_exp = 0.0 # background exp.
                fg_exp = 0.0 # foreground exp.
                # Get sample IDs,
                # avg. expression levels
                # and indices
                sampleIDs = str(
                    tss.qualifiers["sampleIDs"]
                ).split(",")
                avg_exps = str(
                    tss.qualifiers["avg_expression_levels"]
                ).split(",")
                idxs = list(reversed(range(len(sampleIDs))))
                # For each sample ID...
                for i in idxs:
                    if (
                        sampleIDs[i].isdigit() and \
                        isfloat.match(avg_exps[i])
                    ):
                        sampleIDs[i] = int(sampleIDs[i])
                        avg_exps[i] = float(avg_exps[i])
                    else:
                        sampleIDs.pop(i)
                        avg_exps.pop(i)
                # For each sampleID
                for i in range(len(sampleIDs)):
                    # Initialize
                    sampleID = sampleIDs[i]
                    avg_exp = avg_exps[i]
                    # If sample...
                    if sampleID in uid2name:
                        expression.setdefault(
                            uid2name[sampleID],
                            [avg_exp, None]
                        )
                        bg_exp += avg_exp
                        if uid2name[sampleID] in samples:
                            fg_exp += avg_exp
                # If expressed...
                if bg_exp > 0:
                    # For each sample...
                    for sample in expression:
                        expression[sample][1] = \
                            expression[sample][0] * 100.0 / bg_exp
                    # If selectively expressed...
                    if (fg_exp * 100.0) / bg_exp >= percent_exp:
                        tss.qualifiers.setdefault(
                            "expression",
                            expression
                        )
                        tss.score = "%.3f" % \
                            round(
                                (fg_exp * 100.0) / bg_exp, 3
                            )
                        sel_exp_tss.append(tss)

    # Sort TSSs by percentile expression
    sel_exp_tss.sort(
        key=lambda x: float(x.score),
        reverse=True
    )

    return sel_exp_tss

#-------------#
# Main        #
#-------------#

if __name__ == "__main__": main()