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

# Import from GUD module
from GUD import GUDglobals
from GUD.ORM.chrom import Chrom
from GUD.ORM.experiment import Experiment
from GUD.ORM.region import Region
from GUD.ORM.sample import Sample
from GUD.ORM.source import Source

#-------------#
# Functions   #
#-------------#

def parse_args():
    """
    This function parses arguments provided via the command
    line using argparse.
    """

    parser = argparse.ArgumentParser(description="initializes a GUD database for the given genome.")

    parser.add_argument("genome", help="genome assembly")

    # MySQL args
    mysql_group = parser.add_argument_group("mysql arguments")
    mysql_group.add_argument("-d", "--db",
        help="database name (default = given genome assembly)")
    mysql_group.add_argument("-H", "--host", default="localhost",
        help="host name (default = localhost)")
    mysql_group.add_argument("-p", "--passwd",
        help="Password (default = ignore this option)")
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

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()