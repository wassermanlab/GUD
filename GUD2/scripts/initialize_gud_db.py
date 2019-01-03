#!/usr/bin/env python

import re
import argparse
import getpass
from sqlalchemy import create_engine
from sqlalchemy.orm import scoped_session, sessionmaker
from sqlalchemy_utils import create_database, database_exists

# Import from GUD module
from GUD2 import GUDglobals
from GUD2.ORM.chrom import Chrom
from GUD2.ORM.experiment import Experiment
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

    parser = argparse.ArgumentParser(description="initializes a GUD database for the given genome.")

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

    # Initialize GUD database
    initialize_gud_db(args.user, args.host,
        args.port, args.db, args.genome)

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

    if not engine.has_table("chroms"):
        # Initialize
        rows = []
        # Create table
        table = Chrom()
        table.metadata.bind = engine
        table.metadata.create_all(engine)
        # Get UCSC FTP dir/file
        directory, file_name = GUDglobals.get_ucsc_ftp_dir_and_file(
            genome, "chrom_size")
        # For each line...
        for line in GUDglobals.fetch_lines_from_ucsc_ftp_file(
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

    if not engine.has_table("region"):
        # Create table
        table = Region()
        table.metadata.bind = engine
        table.metadata.create_all(engine)

    if not engine.has_table("sample"):
        # Create table
        table = Sample()
        table.metadata.bind = engine
        table.metadata.create_all(engine)

    if not engine.has_table("source"):
        # Create table
        table = Source()
        table.metadata.bind = engine
        table.metadata.create_all(engine)

    if not engine.has_table("experiment"):
        # Create table
        table = Experiment()
        table.metadata.bind = engine
        table.metadata.create_all(engine)

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()