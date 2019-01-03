#!/usr/bin/env python

import re, sys
import argparse
from binning import assign_bin
import getpass
from sqlalchemy import create_engine
from sqlalchemy.orm import scoped_session, sessionmaker
from sqlalchemy_utils import create_database, database_exists

# Import from GUD module
from GUD2 import GUDglobals
from GUD2.ORM.gene import Gene
from GUD2.ORM.region import Region
from GUD2.ORM.source import Source

#-------------#
# Functions   #
#-------------#

def parse_args():
    """
    This function parses arguments provided via the command
    line using argparse.
    """

    parser = argparse.ArgumentParser(description="inserts the UCSC's \"refGene\" table for given genome into GUD.")

    parser.add_argument("genome", help="Genome assembly")

    # MySQL args
    mysql_group = parser.add_argument_group("mysql arguments")
    mysql_group.add_argument("-d", "--db",
        help="Database name (default = given genome assembly)")
    mysql_group.add_argument("-H", "--host", default="localhost",
        help="Host name (default = localhost)")
    mysql_group.add_argument("-P", "--port", default=5506, type=int,
        help="Port number (default = 5506)")

    user = getpass.getuser()
    mysql_group.add_argument("-u", "--user", default=user,
        help="User name (default = current user)")

    args = parser.parse_args()

    # Set default
    if not args.db:
        args.db = args.genome

    return args

def main():

    # Parse arguments
    args = parse_args()

    # Insert genes to GUD database
    insert_genes_to_gud_db(args.user, args.host,
        args.port, args.db, args.genome)

def insert_genes_to_gud_db(user, host, port, db, genome):

    # Initialize
    db_name = "mysql://{}:@{}:{}/{}".format(
        user, host, port, db)
    if not database_exists(db_name):
        create_database(db_name)
    session = scoped_session(sessionmaker())
    engine = create_engine(db_name, echo=False)
    session.remove()
    session.configure(bind=engine, autoflush=False,
        expire_on_commit=False)

    # Get source
    source = Source()
    sou = source.select_by_name(session, "refGene")
    if not sou:
        # Insert source name
        source.name = "refGene"
        session.add(source)
        session.commit()
        sou = source.select_by_name(session, "refGene")

    if not engine.has_table("gene"):
        # Create table
        table = Gene()
        table.metadata.bind = engine
        table.metadata.create_all(engine)
        # Get UCSC FTP dir/file
        directory, file_name = GUDglobals.get_ucsc_ftp_dir_and_file(
            genome, "gene")
        # For each line...
        for line in GUDglobals.fetch_lines_from_ucsc_ftp_file(
            genome, directory, file_name):
            # Split line
            line = line.split("\t")
            # Ignore non-standard chroms, scaffolds, etc.
            m = re.search("^chr(\S+)$", line[2])
            if not m.group(1) in GUDglobals.chroms: continue
            # Get coordinates
            chrom = line[2]
            start = int(line[4])
            end = int(line[5])
            # Get region
            region = Region()
            reg = region.select_unique(
                session, chrom, start, end)
            if not reg:
                # Insert region
                region.bin = assign_bin(start, end)
                region.chrom = chrom
                region.start = start
                region.end = end
                session.add(region)
                session.commit()
                reg = region.select_unique(
                    session, chrom, start, end)
            # Add gene
            gene = Gene()
            if gene.is_unique(session, reg.uid, sou.uid, line[1], line[3]):
                gene.regionID = reg.uid
                gene.sourceID = sou.uid
                gene.name = line[1]
                gene.name2 = line[12]
                gene.strand = line[3]
                gene.cdsStart = line[6]
                gene.cdsEnd = line[7]
                # Using python 3
                if sys.version_info[0] == 3:
                    gene.exonStarts = bytes(line[9], "utf8")
                    gene.exonEnds = bytes(line[10], "utf8")
                else:
                    gene.exonStarts = line[9]
                    gene.exonEnds = line[10]
                session.add(gene)
                session.commit()

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()