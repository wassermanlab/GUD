#!/usr/bin/env python

import argparse
from binning import assign_bin
from ftplib import FTP
import getpass
import gzip
from io import BytesIO
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

#-------------#
# Functions   #
#-------------#

def parse_args():
    """
    This function parses arguments provided via the command
    line using argparse.
    """

    parser = argparse.ArgumentParser(description="this script inserts \"gene\" definitions from the UCSC's \"RefSeq\" table into GUD.")

    parser.add_argument("genome", help="genome assembly")

    # MySQL args
    mysql_group = parser.add_argument_group("mysql arguments")
    mysql_group.add_argument("-d", "--db",
        help="database name (default = given genome assembly)")
    mysql_group.add_argument("-H", "--host", default="localhost",
        help="host name (default = localhost)")
    mysql_group.add_argument("-p", "--passwd",
        help="password (default = ignore this option)")
    mysql_group.add_argument("-P", "--port", default=5506, type=int,
        help="port number (default = 5506)")

    user = getpass.getuser()
    mysql_group.add_argument("-u", "--user", default=user,
        help="user name (default = current user)")

    args = parser.parse_args()

    # Set default
    if not args.db:
        args.db = args.genome
    if not args.passwd:
        args.passwd = ""

    return args

def main():

    # Parse arguments
    args = parse_args()

    # Initialize GUD database: create tables
    # and download data to populate them 
    initialize_gud_db(args.user, args.passwd, args.host,
        args.port, args.db, args.genome)

def initialize_gud_db(user, passwd, host, port, db, genome):

    # Initialize
    db_name = "mysql://{}:{}@{}:{}/{}".format(
        user, passwd, host, port, db)
    if not database_exists(db_name):
        raise ValueError("GUD database does not exist!!!\n\t%s" % db_name)
    session = scoped_session(sessionmaker())
    engine = create_engine(db_name, echo=False)
    session.remove()
    session.configure(bind=engine, autoflush=False,
        expire_on_commit=False)

    table = Gene()
    if not engine.has_table(table.__tablename__):
        # Initialize
        rows = []
        # Create table
        table.__table__.create(bind=engine)
        # Get UCSC FTP file
        directory, file_name = get_ftp_dir_and_file(genome, "gene")
        # Get source
        source = Source()
        sou = source.select_by_name(session, "refGene")
        if not sou: 
            source.name = "refGene"
            session.merge(source)
            session.commit()
            sou = source.select_by_name(session, "refGene")
        # Download data
        for line in fetch_lines_from_ftp_file(
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
            strand = line[3]
            # Get region
            region = Region()
            if region.is_unique(session, chrom, start, end, strand):
                # Insert region
                region.bin = assign_bin(start, end)
                region.chrom = chrom
                region.start = start
                region.end = end
                region.strand = strand
                session.add(region)
                session.commit()
            reg = region.select_unique(session, chrom, start, end, strand)
            # Insert gene
            gene = Gene()
            gene.regionID = reg.uid
            gene.sourceID = sou.uid
            gene.name = line[1]
            gene.name2 = line[12]
            gene.cdsStart = line[6]
            gene.cdsEnd = line[7]
#            gene.strand = line[3]
            gene.exonStarts = line[9]
            gene.exonEnds = line[10]
            session.merge(gene)
            session.commit()

def get_ftp_dir_and_file(genome, data_type):

    # Initialize
    ftp = FTP("hgdownload.soe.ucsc.edu")
    ftp.login()

    # Change into "genome" folder
    try:
        ftp.cwd(os.path.join("goldenPath", genome))
    except:
        raise ValueError("Cannot connect to FTP goldenPath folder: %s" % genome)

    # Fetch bigZips and database files
    if data_type == "gene":
        return "database", "refGene.txt.gz"

def fetch_lines_from_ftp_file(genome, directory, file_name):
    
    # Initialize
    global BIO
    ftp = FTP("hgdownload.soe.ucsc.edu")
    ftp.login()
    BIO = BytesIO()

    # Change into "genome" "directory" folder
    try:
        ftp.cwd(os.path.join("goldenPath", genome, directory))
    except:
        raise ValueError("Cannot connect to FTP goldenPath folder: %s/%s" % (genome, directory))

    # If valid file...
    if file_name in ftp.nlst():
        # Retrieve FTP file
        ftp.retrbinary("RETR %s" % file_name, callback=handle_bytes)
        BIO.seek(0) # Go back to the start
        # If compressed file...
        if file_name.endswith(".gz"):
            f = gzip.GzipFile(fileobj=BIO, mode="rb")
        # ... Else...
        else:
            f = BIO
        # For each line...
        for line in f:
            yield line.decode("UTF-8").strip("\n")

def handle_bytes(bytes):
    BIO.write(bytes)

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()