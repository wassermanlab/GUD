#!/usr/bin/env python

import argparse
import os
import sys

# Import from GUD
from GUD import GUDglobals
from GUD.ORM.chrom import Chrom
from GUD.ORM.gene import Gene
from GUD.ORM.genomic_feature import GenomicFeature
from GUD.ORM.tad import TAD

__doc__ = """
usage: gene2region.py (--gene [STR ...] | --gene-file FILE)
                      [-h] [-d] [--dummy-dir DIR] [-o FILE]
                      [[--sample [STR ...] | --sample-file FILE]]
                      [-d STR] [-H STR] [-p STR] [-P STR] [-u STR]

delimits a region for the given genes based on TADs, nearby genes,
or distance (in kb).

  --gene [STR ...]    gene(s) (e.g. "CD19")
  --gene-file FILE    file containing a list of genes

optional arguments:
  -h, --help          show this help message and exit
  -d, --delimit-by    delimit the gene region based on TADs (i.e.
                      "tad"; default), nearby genes (i.e. "gene"),
                      or distance (i.e. +/- N kb)
  --dummy-dir DIR     dummy directory (default = "/tmp/")
  -o FILE             output file (default = stdout)
  --sample [STR ...]  sample(s) for GUD features (e.g. "B cell")
  --sample-file FILE  file containing a list of samples for GUD
                      features

mysql arguments:
  -d STR, --db STR    database name (default = "%s")
  -H STR, --host STR  host name (default = "%s")
  -p STR, --pwd STR   password (default = ignore this option)
  -P STR, --port STR  port number (default = %s)
  -u STR, --user STR  user name (default = "%s")
""" % \
(
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
    
    # Optional args
    optional_group = parser.add_argument_group(
        "optional arguments"
    )
    optional_group.add_argument(
        "-h", "--help",
        action="store_true"
    )
    optional_group.add_argument(
        "-d", "--delimit-by",
        default="tad"
    )
    optional_group.add_argument(
        "--dummy-dir",
        default="/tmp/"
    )
    optional_group.add_argument("-o")

    # Sample args
    sample_group = parser\
        .add_mutually_exclusive_group()
    sample_group.add_argument(
        "--sample", nargs="*"
    )
    sample_group.add_argument("--sample-file")

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

    if "-h" in sys.argv or "--help" in sys.argv:
        print(__doc__)
        exit(0)

    # Mandatory arguments
    gene_group = parser\
        .add_mutually_exclusive_group(
            required=True
        )
    gene_group.add_argument(
        "--gene", nargs="*"
    )
    gene_group.add_argument("--gene-file")

    args = parser.parse_args()

    if not args.delimit_by.isdigit():
        if args.delimit_by not in [
            "tad", "gene"
        ]:
            print(__doc__)
            exit(0)

    return args

def main():

    # Parse arguments
    args = parse_args()

    # Initialize
    dummy_file = os.path.join(
        os.path.abspath(args.dummy_dir),
        "%s.%s.bed" % (
            os.path.basename(__file__),
            os.getpid()
        )
    )

    # Fetch genes
    genes = []
    if args.gene:
        genes = args.gene
    else:
        gene_file = os.path.abspath(
            args.gene_file)
        for g in GUDglobals.parse_file(
            gene_file):
            genes.append(g)

    # Fetch samples
    samples = []
    if args.sample:
        samples = args.sample
    elif args.sample_file:
        sample_file = os.path.abspath(
            args.sample_file)
        for s in GUDglobals.parse_file(
            sample_file
        ):
            samples.append(s)

    # Establish SQLalchemy session with GUD
    session = GUDglobals.establish_GUD_session(
        args.user,
        args.pwd,
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
            args.delimit_by
        )
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
        for line in GUDglobals.parse_file(
            dummy_file
        ):
            GUDglobals.write(None, line)

    # Delete dummy file
    os.remove(dummy_file)
    
def get_gene_region(session, gene, samples=[],
    delimit_by="tad"):
    """
    Delimits a region for the given gene based
    on TAD boundaries, nearby genes, or +/- N
    kb.
    """

    # Initialize
    chrom = None

    # Get chromosome sizes
    chrom_sizes = Chrom.chrom_sizes(session)

    # Get all genes with the given name
    genes = Gene.select_by_name(session,
        gene, as_genomic_feature=True)
    # If genes are not valid...
    if not genes:
        raise ValueError(
            "Gene \"%s\" is not valid!!!" % gene
        )

    # For each gene...
    for g in genes:

        # Ignore non-standard chroms,
        # scaffolds, etc.
        if g.chrom in chrom_sizes:

            # If first gene...
            if not chrom:
                chrom = g.chrom
                gene_start = chrom_sizes[chrom]
                gene_end = 0
                region_start = 0
                region_end = chrom_sizes[chrom]

            # If gene is on a different chrom...
            if g.chrom != chrom:
                raise ValueError(
                    "Gene \"%s\" is mapped to different chromosomes!!!" % gene
                )

            # If start position is upstream
            # from previous...
            if g.start < gene_start:
                gene_start = g.start

            # If end position is downstream
            # from previous...
            if g.end > gene_end:
                gene_end = g.end

    # If chrom...
    if chrom:

        # If delimit by distance...
        if delimit_by.isdigit():
            region_start = gene_start -\
                int(delimit_by) * 1000
            region_end = gene_end +\
                int(delimit_by) * 1000

        # Instead, if delimit by TADs...
        elif delimit_by == "tad":
            region_start, region_end =\
                get_region_coordinates_by_tad(
                    session,
                    chrom,
                    gene_start,
                    gene_end,
                    region_start,
                    region_end,
                    samples
                )

        # Instead, if delimit by genes...
        elif delimit_by == "gene":
            region_start, region_end =\
                get_region_coordinates_by_gene(
                    session,
                    gene,
                    chrom,
                    gene_start,
                    gene_end,
                    region_start,
                    region_end
                )

        # Get genomic feature of region
        gene_region = GenomicFeature(
            chrom,
            region_start,
            region_end,
            feat_type = "Region",
            feat_id = "%s" % gene
        )

        return gene_region

    raise ValueError(
        "Gene \"%s\" is not in a valid chromosome!!!" % gene
    )

#def get_region_coordinates_by_tad(session,
#    chrom, gene_start, gene_end, region_start,
#    region_end, samples=[]):
#    """
#    Delimits a region for the given gene based
#    on TAD boundaries.
#    """
#
#    # Get TADs encompassing the gene
#    tads = Tad.select_by_location(
#        session,
#        chrom,
#        gene_start,
#        gene_end,
#        samples,
#        as_genomic_feature = True
#    )
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
    """
    Delimits a region for the given gene based
    on nearby genes.
    """

    # Get genes within chromosome
    genes = Gene.select_by_location(
        session,
        chrom,
        region_start,
        region_end,
        as_genomic_feature = True
    )

    # For each gene...
    for g in genes:
        # Skip the gene
        if g.id == gene: continue
        # Get the closest upstream gene end
        # to the gene's start
        if g.end <= gene_start\
        and g.end > region_start:
            region_start = g.end
        # Get the closest downstream gene
        # start to the gene's end
        if g.start >= gene_end\
        and g.start < region_end:
            region_end = g.start

    return region_start, region_end

#-------------#
# Main        #
#-------------#

if __name__ == "__main__": main()