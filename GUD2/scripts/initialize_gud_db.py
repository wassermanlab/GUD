#!/usr/bin/env python

import os, sys, re
import argparse
from binning import assign_bin
from datetime import date
from ftplib import FTP
import getpass
import gzip
from io import BytesIO
from sqlalchemy import create_engine
from sqlalchemy.orm import scoped_session, sessionmaker
from sqlalchemy_utils import create_database, database_exists

from GUD2 import GUDglobals
from GUD2.ORM.chroms import Chroms
from GUD2.ORM.region import Region
from GUD2.ORM.sample import Sample
from GUD2.ORM.source import Source
from GUD2.ORM.experiment import Experiment

#-------------#
# Functions   #
#-------------#

def parse_args():
    """
    This function parses arguments provided via the command
    line using argparse.
    """
    parser = argparse.ArgumentParser(description="this script initializes a GUD database for the given genome.")
    parser.add_argument("genome", help="Genome assembly")

    # MySQL args
    mysql_group = parser.add_argument_group("mysql arguments")
    mysql_group.add_argument("-d", "--db",
        help="Database name (default = input genome assembly)")
    mysql_group.add_argument("-H", "--host", default="localhost",
        help="Host name (default = localhost)")
    mysql_group.add_argument("-P", "--port", default=5506, type=int,
        help="Port number (default = 5506)")
    mysql_group.add_argument("-u", "--user", default=getpass.getuser(),
        help="User name (default = current user)")

    args = parser.parse_args()

    # Set default
    if not args.db:
        args.db = args.genome

    return args

def initialize_gud_db(user, host, port, db, genome):

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
    today = str(date.today())

    # Create chrom sizes table
    if not engine.has_table("chroms"):
        # Initialize
        rows = []
        table = Chroms()
        table.metadata.bind = engine
        # Create table
        table.metadata.create_all(engine)
        # Get UCSC FTP file
        directory, file_name = get_ftp_dir_and_file(genome, "chrom_size")
        # Download data
        for line in fetch_lines_from_ftp_file(
            genome, directory, file_name):
            # Split line
            line = line.split("\t")
            # Ignore non-standard chroms, scaffolds, etc.
            print line
            m = re.search("^chr(\S+)$", line[0])
	    if not m.group(1) in GUDglobals.chroms: continue
            # Add row
            rows.append(
                {
                    "chrom": line[0],
                    "size": line[1]
                }
            )
        # Insert rows to table
        engine.execute(table.__table__.insert(), rows)
    if not engine.has_table("region"):
        # initialize table
        table = Region()
        table.metadata.bind = engine
        table.metadata.create_all(engine)
    if not engine.has_table("sample"):
        # initialize table
        table = Sample()
        table.metadata.bind = engine
        table.metadata.create_all(engine)
    if not engine.has_table("source"):
        # initialize table
        table = Source()
        table.metadata.bind = engine
        table.metadata.create_all(engine)
    if not engine.has_table("experiment"):
        # initialize table
        table = Experiment()
        table.metadata.bind = engine
        table.metadata.create_all(engine)

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
    if data_type == "chrom_size":
        return "bigZips", "%s.chrom.sizes" % genome

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

    # Parse arguments
    args = parse_args()

    # Initialize GUD database: create tables
    # and download data to populate them 
    initialize_gud_db(args.user, args.host,
        args.port, args.db, args.genome)
