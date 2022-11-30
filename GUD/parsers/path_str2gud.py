#!/usr/bin/env python
# mysql -u ontarget_w --database tamar_test
import os
import sys
import re
import argparse
from binning import assign_bin
import getpass
from sqlalchemy import create_engine
from sqlalchemy.orm import mapper, scoped_session, sessionmaker
from sqlalchemy_utils import database_exists
import warnings

# Import from GUD module
from GUD import GUDglobals
from GUD.ORM.short_tandem_repeat import ShortTandemRepeat
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

    parser = argparse.ArgumentParser(
        description="this script inserts short tandem repeat location information")

    parser.add_argument("bed_file", help="gangSTR bed file")

    # Optional args
    parser.add_argument("--source", default="gangSTR", help="Source name")

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


def insert_str_to_gud_db(user, host, port, db, bed_file, source_name):

    # Initialize
    metadata = {}
    db_name = "mysql://{}:@{}:{}/{}".format(
        user, host, port, db)
    if not database_exists(db_name):
        raise ValueError("GUD db does not exist: %s" % db_name)
    session = scoped_session(sessionmaker())
    engine = create_engine(db_name, echo=False)
    session.remove()
    session.configure(bind=engine, autoflush=True,
                      expire_on_commit=False)

    # Initialize table
    if not engine.has_table("short_tandem_repeats"):
        raise ValueError("No str table!")
    if not engine.has_table("regions"):
        raise ValueError("No regions table!")
    if not engine.has_table("sources"):
        raise ValueError("No sources table!")

    # source entry
    source = Source()
    sou = source.select_by_name(session, source_name)
    if not sou:
            source.name = source_name
            session.merge(source)
            session.commit()
            sou = source.select_by_name(session, source_name)

    # parse table
    with open(bed_file) as f:
        for line in f:
            if not line.startswith("#"):
                split_line = line.split("\t")
                chrom = str(split_line[0])
                start = int(split_line[1])
                end = int(split_line[2])
                motif = str(split_line[3])
                pathogenicity = int(split_line[4].rstrip())

                # region entry
                region = Region()
                reg = region.select_unique(
                    session, chrom, start, end)
                if not reg:
                    region.bin = assign_bin(start, end)
                    region.chrom = chrom
                    region.start = start
                    region.end = end
                    session.merge(region)
                    session.commit()
                    reg = region.select_unique(
                        session, chrom, start, end)

                # str entry
                STR = ShortTandemRepeat()
                if STR.is_unique(session, reg.uid, sou.uid, pathogenicity):
                    STR.motif = motif
                    STR.pathogenicity = pathogenicity
                    STR.regionID = reg.uid
                    STR.sourceID = sou.uid
                    session.merge(STR)
                    session.commit()

#-------------#
# Main        #
#-------------#


if __name__ == "__main__":

    # Parse arguments
    args = parse_args()

    # Insert ENCODE data to GUD database
    insert_str_to_gud_db(args.user, args.host, args.port, args.db,
                         args.bed_file, args.source)
