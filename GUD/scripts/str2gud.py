#!/usr/bin/env python

import os, sys, re
import argparse
from binning import assign_bin
from datetime import date
import getpass
from sqlalchemy import create_engine
from sqlalchemy.orm import mapper, scoped_session, sessionmaker
from sqlalchemy_utils import database_exists
import warnings

# Import from GUD module
from GUD import GUDglobals
from GUD.ORM.short_tandem_repeat import ShortTandemRepeat

#-------------#
# Classes     #
#-------------#

class Model(object): pass

#-------------#
# Functions   #
#-------------#

def parse_args():
    """
    This function parses arguments provided via the command
    line using argparse.
    """

    parser = argparse.ArgumentParser(description="this script inserts short tandem repeat location information")

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
    session.configure(bind=engine, autoflush=False,
        expire_on_commit=False)
    today = str(date.today())

    # Initialize table
    table = ShortTandemRepeat()
    table.metadata.bind = engine
    # if not engine.has_table("short_tandem_repeat"):
    try:
        table.metadata.create_all(engine)
    except:
        raise ValueError("Cannot create \"short_tandem_repeat\" table!")
    mapper(Model, table.__table__)

    # parse table 
    with open(bed_file) as f:
        for line in f: 
            split_line = line.split("\t")
            chrom = str(split_line[0])
            start = int(split_line[1])
            end = int(split_line[2])
            length = int(split_line[3])
            motif = str(split_line[4])
            pathogenicity = 0

            model = Model()
            model.bin             = assign_bin(start, end)
            model.chrom           = chrom
            model.start           = start
            model.end             = end
            model.length          = length
            model.motif           = motif
            model.pathogenicity   = pathogenicity
            model.date = today

            session.merge(model)
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
