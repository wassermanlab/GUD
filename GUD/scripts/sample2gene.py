#!/usr/bin/env python

import argparse
import os
import shutil

# Import from GUD
from GUD import GUDglobals
from GUD.ORM.chrom import Chrom
from GUD.ORM.expression import Expression
from GUD.ORM.gene import Gene
from GUD.ORM.genomic_feature import GenomicFeature
from GUD.ORM.sample import Sample
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
        help="require expression in all sample(s) (default = False)")
    parser.add_argument("--max-genes", metavar="", type=int,
        default=GUDglobals.max_genes,
        help="max. number of genes to return (default = %s)" % GUDglobals.max_genes)
    parser.add_argument("--min-exp", metavar="", type=float,
        default=GUDglobals.min_exp,
        help="min. expression in input sample(s) (in TPM; default = %s)" % GUDglobals.min_exp)
    parser.add_argument("--min-percent-exp", metavar="", type=float,
        default=GUDglobals.min_percent_exp,
        help="min. percentage of expression in input sample(s) (default = %s)" % GUDglobals.min_percent_exp)

    # MySQL args
    mysql_group = parser.add_argument_group("mysql arguments")
    mysql_group.add_argument("-d", "--db", default=GUDglobals.db_name,
        help="database name (default = \"%s\")" % GUDglobals.db_name)
    mysql_group.add_argument("-H", "--host", default=GUDglobals.db_host,
        help="host name (default = \"%s\")" % GUDglobals.db_host)
    mysql_group.add_argument("-p", "--pwd", metavar="PASS",
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
        for sample in GUDglobals.parse_file(sample_file):
            samples.append(sample)
    else:
        raise ValueError("No sample(s) was provided!")

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
        args.min_exp,
        args.min_percent_exp,
        args.all
    )

    # If differentially expressed TSSs...
    if diff_exp_tss:

        # For each TSS...
        for tss in diff_exp_tss[:args.max_genes]:
            # Write
            GUDglobals.write(dummy_file, tss)

        # If output file...
        if args.out_file:
            shutil.copy(
                dummy_file,
                os.path.abspath(args.out_file)
            )
        # ... Else, print on stdout...
        else:
            for line in GUDglobals.parse_file(dummy_file):
                GUDglobals.write(None, line)

        # Delete dummy file
        os.remove(dummy_file)

    else:
        raise ValueError("No differentially expressed gene(s) found!!!")

def get_differentially_expressed_tss(session,
    samples=[], min_tpm=100, min_percent_exp=25,
    exp_in_all_samples=False):
    """
    Identifies TSSs differentially expressed in
    samples.
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
        session, samples, min_tpm):
        # Get index
        i = samples.index(tss.Sample.name)
        sample2tss[i].add(tss.Expression.tssID)
        tssIDs.add(tss.Expression.tssID)
    # If required expression in all samples...
    if exp_in_all_samples:
        tssIDs = list(set.intersection(*sample2tss))
    # ... Else...
    else:
        tssIDs = list(tssIDs)

    # If there are TSS IDs...
    if tssIDs:
        # For each TSS...
        for tss in TSS.select_by_uids(session,
            tssIDs, as_genomic_feature=True):
            # Initialize
            expression = {}
            background_exp = 0.0
            foreground_exp = 0.0
            # For each sampleID
            for i in range(len(tss.qualifiers["sampleIDs"])):
                if tss.qualifiers["sampleIDs"][i] in uid2name:
                    expression.setdefault(
                        uid2name[tss.qualifiers["sampleIDs"][i]],
                        tss.qualifiers["avg_expression_levels"][i]
                    )
                    background_exp += tss.qualifiers["avg_expression_levels"][i]
                    if uid2name[tss.qualifiers["sampleIDs"][i]] in samples:
                        foreground_exp += tss.qualifiers["avg_expression_levels"][i]
            # If expressed...
            if background_exp > 0:
                # Get percent exp. in samples
                percent_exp = (foreground_exp * 100.0) / background_exp
                # If differentially expressed...
                if percent_exp >= min_percent_exp:
                    tss.qualifiers.setdefault("expression", expression)
                    tss.qualifiers.setdefault("percent_exp", percent_exp)
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