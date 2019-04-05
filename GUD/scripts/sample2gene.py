#!/usr/bin/env python

import argparse
import os
import shutil
import sys

# Import from GUD
from GUD import GUDglobals
from GUD.ORM.chrom import Chrom
from GUD.ORM.expression import Expression
from GUD.ORM.gene import Gene
from GUD.ORM.genomic_feature import GenomicFeature
from GUD.ORM.sample import Sample
from GUD.ORM.tss import TSS

usage_msg = """
usage: sample2gene.py (--sample [STR ...] | --sample-file FILE)
                      [-h] [--dummy-dir DIR] [-o FILE]
                      [-a] [--percent FLT] [--tpm FLT] [--tss INT]
                      [-d STR] [-H STR] [-p STR] [-P STR] [-u STR]
"""

help_msg = """%s

identifies one or more genes differentially expressed in samples.

  --sample [STR ...]  sample(s) (e.g. "B cell")
  --sample-file FILE  file containing a list of samples

optional arguments:
  -h, --help          show this help message and exit
  --dummy-dir DIR     dummy directory (default = "/tmp/")
  -o FILE             output file (default = stdout)

expression arguments:
  -a, --all           expression in all samples (default = False)
  --percent FLT       min. percentage of expression for TSS in
                      input samples (default = ignore this option)
  --tpm FLT           min. expression levels (in TPM) for TSS in
                      input samples (default = %s)
  --tss INT           max. number of TSSs to return (default = %s)

mysql arguments:
  -d STR, --db STR    database name (default = "%s")
  -H STR, --host STR  host name (default = "%s")
  -p STR, --pwd STR   password (default = ignore this option)
  -P STR, --port STR  port number (default = %s)
  -u STR, --user STR  user name (default = "%s")
""" % \
(
    usage_msg,
    GUDglobals.min_tpm_exp,
    GUDglobals.max_num_tss,
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
        default=0.0
    )
    exp_group.add_argument(
        "--tpm",
        default=GUDglobals.min_tpm_exp
    )
    exp_group.add_argument(
        "--tss",
        default=GUDglobals.max_num_tss
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

    # Get TSSs
    diff_exp_tss = get_differentially_expressed_tss(
        session,
        samples,
        args.tpm,
        args.percent,
        args.all
    )

    # If differentially expressed TSSs...
    if diff_exp_tss:

        # For each TSS...
        for tss in diff_exp_tss[:args.tss]:
            # Write
            GUDglobals.write(dummy_file, tss)

        # If output file...
        if args.o:
            shutil.copy(
                dummy_file,
                os.path.abspath(args.o)
            )
        # ... Else, print on stdout...
        else:
            for line in GUDglobals.parse_file(dummy_file):
                GUDglobals.write(None, line)

        # Delete dummy file
        os.remove(dummy_file)

    else:
        raise ValueError(
            "No differentially expressed genes found!!!"
        )

def get_differentially_expressed_tss(session,
    samples=[], tpm_exp=100, percent_exp=0.0,
    exp_in_all_samples=False):
    """
    Identifies TSSs from genes differentially
    expressed in samples.
    """

    # Initialize
    name2uid = {}
    uid2name = {}
    tssIDs = set()
    sample2tss = []
    diff_exp_tss = []

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
        session, samples, tpm_exp):
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
                # For each sampleID
                for i in range(
                    len(tss.qualifiers["sampleIDs"])
                ):
                    # Initialize
                    sampleID = tss\
                        .qualifiers["sampleIDs"][i]
                    avg_exp = tss\
                        .qualifiers["avg_expression_levels"][i]
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
                    # If differentially expressed...
                    if (fg_exp * 100.0) / bg_exp >= percent_exp:
                        tss.qualifiers.setdefault(
                            "expression",
                            expression
                        )
                        tss.qualifiers.setdefault(
                            "percent_exp",
                            (fg_exp * 100.0) / bg_exp
                        )
                        diff_exp_tss.append(tss)

    # Sort TSSs by percent exp.
    diff_exp_tss.sort(
        key=lambda x: x.qualifiers["percent_exp"],
        reverse=True
    )

    return diff_exp_tss

#-------------#
# Main        #
#-------------#

if __name__ == "__main__": main()