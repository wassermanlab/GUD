#!/usr/bin/env python

import argparse
from binning import assign_bin
import getpass
import os
import re
from sqlalchemy import create_engine
from sqlalchemy.orm import (
    scoped_session,
    sessionmaker
)
from sqlalchemy_utils import database_exists

# Import from GUD module
from GUD import GUDglobals
from GUD.ORM.gene import Gene
from GUD.ORM.region import Region
from GUD.ORM.source import Source
from .initialize import (
    initialize_gud_db,
    get_ftp_dir_and_file,
    fetch_lines_from_ftp_file
)

usage_msg = """
usage: genes2gud.py --genome STR [-h]
                    [-d STR] [-H STR] [-p STR] [-P STR] [-u STR]
"""

help_msg = """%s

inserts gene definitions from the UCSC's RefSeq table into GUD.

  --genome STR        genome assembly

optional arguments:
  -h, --help          show this help message and exit

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

    # Mandatory arguments
    parser.add_argument("--genome")

    # Optional args
    optional_group = parser.add_argument_group(
        "optional arguments"
    )
    optional_group.add_argument(
        "-h", "--help",
        action="store_true"
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

    # Print help
    if args.help:
        print(help_msg)
        exit(0)

    # Check mandatory arguments
    if (not args.genome):
        print(": "\
            .join(
                [
                    "%s\ngenes2gud.py" % usage_msg,
                    "error",
                    "argument \"--genome\" is required\n"
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

    # Insert RefGene data
    refgene_to_gud(
        args.user,
        args.pwd,
        args.host,
        args.port,
        args.db,
        args.genome
    )

def refgene_to_gud(user, pwd, host, port, db,
    genome):

    # Initialize
    features = []
    db_name = "mysql://{}:{}@{}:{}/{}".format(
        user, pwd, host, port, db)
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
        echo=False
    )
    session.remove()
    session.configure(
        bind=engine,
        autoflush=False,
        expire_on_commit=False
    )

    table = Gene()
    if not engine.has_table(
        table.__tablename__
    ):
        # Initialize
        rows = []
        # Create table
        table.__table__.create(bind=engine)
        # Get UCSC FTP file
        directory, file_name =\
            get_ftp_dir_and_file(
                genome,
                "gene"
            )
        # Get source
        source = Source()
        sou = source.select_by_name(
            session,
            "refGene"
        )
        if not sou: 
            source.name = "refGene"
            session.merge(source)
            session.commit()
            sou = source.select_by_name(
                session,
                "refGene"
            )
        # Download data
        for line in fetch_lines_from_ftp_file(
            genome,
            directory,
            file_name
        ):
            # Split line
            line = line.split("\t")
            # Ignore non-standard chroms,
            # scaffolds, etc.
            m = re.search("^chr(\S+)$", line[2])
            if not m.group(1) in GUDglobals.chroms:
                continue
            # Get coordinates
            chrom = line[2]
            start = int(line[4])
            end = int(line[5])
            strand = line[3]
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
            # Insert gene
            gene = Gene()
            gene.regionID = reg.uid
            gene.sourceID = sou.uid
            gene.name = line[1]
            gene.name2 = line[12]
            gene.cdsStart = line[6]
            gene.cdsEnd = line[7]
            gene.exonStarts = line[9]
            gene.exonEnds = line[10]
            session.add(gene)
            session.commit()

#-------------#
# Main        #
#-------------#

if __name__ == "__main__": main()