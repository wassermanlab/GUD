#!/usr/bin/env python

import argparse
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
from sqlalchemy_utils import (
    create_database,
    database_exists
)

from GUD2 import GUDglobals
from GUD2.ORM.chrom import Chrom
from GUD2.ORM.experiment import Experiment
from GUD2.ORM.expression import Expression
from GUD2.ORM.region import Region
from GUD2.ORM.sample import Sample
from GUD2.ORM.source import Source

#-------------#
# Functions   #
#-------------#

def parse_args():
    """
    This function parses arguments provided via the command
    line using argparse.
    """

    parser = argparse.ArgumentParser(description="this script initializes a GUD database for the given genome.")

    parser.add_argument("genome", help="genome assembly")

    # MySQL args
    mysql_group = parser.add_argument_group("mysql arguments")
    mysql_group.add_argument("-d", "--db",
        help="database name (default = given genome assembly)")
    mysql_group.add_argument("-H", "--host", default="localhost",
        help="host name (default = localhost)")
    mysql_group.add_argument("-p", "--passwd",
        help="Password (default = do not use)")
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

    # Initialize GUD database
    initialize_gud_db(args.user, args.passwd, args.host,
        args.port, args.db, args.genome)

def initialize_gud_db(user, passwd, host, port, db, genome):

    # Initialize
    db_name = "mysql://{}:{}@{}:{}/{}".format(
        user, passwd, host, port, db)
    if not database_exists(db_name):
        create_database(db_name)
    session = scoped_session(sessionmaker())
    engine = create_engine(db_name, echo=False)
    session.remove()
    session.configure(bind=engine, autoflush=False,
        expire_on_commit=False)

    table = Chrom()
    if not engine.has_table(table.__tablename__):
        # Intialize
        rows = []
        # Create table
        table.__table__.create(bind=engine)
        # Get UCSC FTP file
        directory, file_name = get_ftp_dir_and_file(genome, "chrom_size")
        # Download data
        for line in fetch_lines_from_ftp_file(
            genome, directory, file_name):
            # Split line
            line = line.split("\t")
            # Ignore non-standard chroms, scaffolds, etc.
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

    table = Experiment()
    if not engine.has_table(table.__tablename__):
        # Create table
        table.__table__.create(bind=engine)

    table = Expression()
    if not engine.has_table(table.__tablename__):
        # Create table
        table.__table__.create(bind=engine)

    table = Region()
    if not engine.has_table(table.__tablename__):
        # Create table
        table.__table__.create(bind=engine)

    table = Sample()
    if not engine.has_table(table.__tablename__):
        # Create table
        table.__table__.create(bind=engine)

    table = Source()
    if not engine.has_table(table.__tablename__):
        # Create table
        table.__table__.create(bind=engine)

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