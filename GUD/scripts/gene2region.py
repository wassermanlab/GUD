#!/usr/bin/env python

import argparse
import os

# Import from GUD
from GUD import GUDglobals
from GUD.ORM.chrom import Chrom
from GUD.ORM.gene import Gene
from GUD.ORM.gud_feature import GUDFeature

#-------------#
# Functions   #
#-------------#

def parse_args():
    """
    This function parses arguments provided via command
    line and returns an {argparse} object.
    """

    parser = argparse.ArgumentParser(description="delimits a gene region.")

    parser.add_argument("--dummy-dir", default="/tmp/",
        help="dummy directory (default = /tmp/)")

    gene_group = parser.add_mutually_exclusive_group(required=True)
    gene_group.add_argument("--gene", nargs="*",
        help="gene symbol(s) (e.g. \"TH\")")
    gene_group.add_argument("--gene-file",
        help="file containing a list of gene symbols")

    parser.add_argument("-l", "--limit-by", default="tad",
        help="limit gene region by: 1) TAD boundaries (i.e. \"tad\"; default); 2) boundaries of nearby genes (i.e. \"gene\"); or 3) +/- N kb around that gene body (e.g. use 1000 for 1 Mb)")

    sample_group = parser.add_mutually_exclusive_group()
    sample_group.add_argument("--sample", default=[], nargs="*",
        help="sample(s) (restrict to features in sample(s); e.g. \"brain\")")
    sample_group.add_argument("--sample-file",
        help="file containing a list of samples")

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
        "%s.%s.bed" % (os.path.basename(__file__),
        os.getpid()))

    # Fetch genes
    genes = []
    if args.gene:
        genes = args.gene
    else:
        gene_file = os.path.abspath(args.gene_file)
        for g in GUDglobals.parse_file(gene_file):
            genes.append(g)

    # Fetch samples
    samples = []
    if args.sample:
        samples = args.sample
    elif args.sample_file:
        sample_file = os.path.abspath(args.sample_file)
        for s in GUDglobals.parse_file(sample_file):
            samples.append(s)

    # Establish SQLalchemy session with GUD
    session = GUDglobals.establish_GUD_session(
        args.user,
        args.passwd,
        args.host,
        args.port,
        args.db
    )

    # For each gene...
    for gene in genes:
        # Get that gene's region
        region = get_gene_region(
            session,
            gene,
            samples,
            args.limit_by)
        # Write
        GUDglobals.write(dummy_file, region)

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
    
def get_gene_region(
    session, gene, sample=[], limit_by="tad"):

    # Initialize
    chrom = None

    # Get chromosome sizes
    chrom_sizes = Chrom.chrom_sizes(session)

    # Fetch all genes with the given gene name
    genes = Gene.select_by_name(
        session, gene, as_gud_feature=True)
    # If genes are not valid...
    if not genes:
        raise ValueError("Gene name is not a valid: %s" % gene)

    # For each gene...
    for g in genes:
        # Ignore non-standard chroms, scaffolds, etc.
        if g.chrom in chrom_sizes:
            # If first gene...
            if not chrom:
                chrom = g.chrom
                gene_start = chrom_sizes[chrom]
                gene_end = 0
                region_start = 0
                region_end = chrom_sizes[chrom]
            # If gene's chromosome is different from previous genes...
            if g.chrom != chrom:
                warnings.warn("\"{}\" chrom (i.e. \"{}\") is different from that of previous genes (i.e. \"{}\")!".format(
                    g.id, g.chrom, chrom))
                return None
            # If start position is upstream from previous start...
            if g.start < gene_start:
                gene_start = g.start
            # If end position is downstream from previous end...
            if g.end > gene_end:
                gene_end = g.end

    # If chrom...
    if chrom:
        # If limit by kb...
        if limit_by.isdigit():
            region_start = \
                gene_start - int(limit_by) * 1000
            region_end = \
                gene_end + int(limit_by) * 1000
        # ... Instead, if limit by TAD...
        elif limit_by == "tad":
            region_start, region_end = \
                get_region_coordinates_by_tad(
                    session,
                    chrom,
                    gene_start,
                    gene_end,
                    region_start,
                    region_end,
                    sample
                )
        # ... Instead, if limit by gene...
        elif limit_by == "gene":
            region_start, region_end = \
                get_region_coordinates_by_gene(
                    session,
                    gene,
                    chrom,
                    gene_start,
                    gene_end,
                    region_start,
                    region_end
                )
        # Define region
        gene_region = GUDFeature(
            chrom,
            region_start,
            region_end,
            feat_type = "GeneRegion",
            feat_id = "%s" % gene
        )

        return gene_region

    warnings.warn("\"{}\" is not on a valid chrom!".format(gene))

    return None

#def get_region_coordinates_by_tad(session, chrom, gene_start, gene_end,
#    region_start, region_end, sample=[]):
#
#    # Get TADs encompassing the gene
#    tads = GUDglobals.Tad.select_encompasing_range(session, chrom, gene_start,
#        gene_end, sample=sample)
#    # If TADs could not be found...
#    if not tads:
#        # Get TADs overlapping the gene
#        tads = GUDglobals.Tad.select_overlapping_range(session, chrom, gene_start,
#            gene_end, sample=sample)
#    # If TADs could not be found...
#    if not tads and sample:
#        warnings.warn("\nCell-/tissue-specific TADs could not be found!\n\tUsing TADs from other cells/tissues instead...\n")
#        # Get any TADs encompassing the gene
#        tads = GUDglobals.Tad.select_encompasing_range(session, chrom,
#            gene_start, gene_end)
#        # If TADs could not be found...
#        if not tads:
#            # Get any TADs overlapping the gene
#            tads = GUDglobals.Tad.select_overlapping_range(session, chrom,
#                gene_start, gene_end)
#
#    # For each TAD... #
#    for t in tads:
#        # Get upstream start closest to the gene's start 
#        if t.start <= gene_start and t.start > region_start:
#            region_start = t.start
#        # Get downstream end closest to the gene's end 
#        if t.end >= gene_end and t.end < region_end:
#            region_end = t.end
#
#    return region_start, region_end

def get_region_coordinates_by_gene(session,
    gene, chrom, gene_start, gene_end,
    region_start, region_end):

    # Get genes within chromosome
    genes = Gene.select_by_location(
        session,
        chrom,
        region_start,
        region_end,
        as_gud_feature = True
    )

    # For each gene...
    for g in genes:
        # Skip the gene
        if g.id == gene: continue
        # Get the closest upstream gene end
        # to the gene's start 
        if g.end <= gene_start and g.end > region_start:
            region_start = g.end
        # Get the closest downstream gene start to the gene's end 
        if g.start >= gene_end and g.start < region_end:
            region_end = g.start

    return region_start, region_end

#-------------#
# Main        #
#-------------#

if __name__ == "__main__": main()