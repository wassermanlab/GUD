#!/usr/bin/env python

import re
import argparse
from binning import assign_bin
import getpass
from sqlalchemy import create_engine
from sqlalchemy.orm import scoped_session, sessionmaker
from sqlalchemy_utils import create_database, database_exists

# Import from GUD module
from GUD2 import GUDglobals
from GUD2.ORM.conservation import Conservation
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

    parser = argparse.ArgumentParser(description="inserts the UCSC's \"multizNway\" table for given genome into GUD.")

    parser.add_argument("genome", help="genome assembly")

    # MySQL args
    mysql_group = parser.add_argument_group("mysql arguments")
    mysql_group.add_argument("-d", "--db",
        help="database name (default = given genome assembly)")
    mysql_group.add_argument("-H", "--host", default="localhost",
        help="host name (default = localhost)")
    mysql_group.add_argument("-P", "--port", default=5506, type=int,
        help="port number (default = 5506)")

    user = getpass.getuser()
    mysql_group.add_argument("-u", "--user", default=user,
        help="user name (default = current user)")

    args = parser.parse_args()

    # Set default
    if not args.db:
        args.db = args.genome

    return args

def main():

    # Parse arguments
    args = parse_args()

    # Insert conservation to GUD database
    insert_conservation_to_gud_db(args.user, args.host,
        args.port, args.db, args.genome)

def insert_conservation_to_gud_db(user, host, port,
    db, genome):

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

    if not engine.has_table("conservation"):
        # Initialize
        rows = []
        # Create table
        table = Conservation()
        table.metadata.bind = engine
        table.metadata.create_all(engine)
        # Get UCSC FTP dir/file
        directory, file_name = GUDglobals.get_ucsc_ftp_dir_and_file(
            genome, "conservation")
        # Get source
        source = Source()
        m = re.search("^(.+).txt.gz$", file_name)
        source_name = m.group(1)
        sou = source.select_by_name(session, source_name)
        if not sou:
            # Insert source name
            source.name = source_name
            session.add(source)
            session.commit()
            sou = source.select_by_name(session, source_name)
        # For each line...
        for line in GUDglobals.fetch_lines_from_ucsc_ftp_file(
            genome, directory, file_name):
            # Split line
            line = line.split("\t")
            # Ignore non-standard chroms, scaffolds, etc.
            m = re.search("^chr(\S+)$", line[1])
            if not m.group(1) in GUDglobals.chroms: continue
            # Get coordinates
            chrom = line[1]
            start = int(line[2])
            end = int(line[3])
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
            # Insert conserved regions in bulks of 100,000
            conservation = Conservation()
            if conservation.is_unique(session, reg.uid, sou.uid):
                conservation.score = line[6]
                conservation.regionID = reg.uid
                conservation.sourceID = sou.uid
                rows.append(conservation)
            if len(rows) == 100000:
                session.add_all(rows)
                session.commit()
                # Initialize
                rows = []
        # Insert remaining conserved regions
        session.add_all(rows)
        session.commit()

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()