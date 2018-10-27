#!/usr/bin/env python

import os, sys, re
import argparse
from binning import assign_bin
from datetime import date
import getpass
from sqlalchemy import create_engine
from sqlalchemy.orm import mapper, scoped_session, sessionmaker
from sqlalchemy_utils import database_exists

# Import from GUD module
from GUD import GUDglobals
from GUD.ORM.tf_binding import TfBinding

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

    parser = argparse.ArgumentParser(description="this script inserts tf data from ReMap into GUD. argument \"directory\" refers to where \"*_TF_archive_all_macs2_*.tar.gz\" was uncompressed.\"")

    parser.add_argument("directory", help="Downloads directory")

    # Optional args
    parser.add_argument("--source", default="ReMap", help="Source name (e.g. \"PMID:29126285\"; default = \"ReMap\")")

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

def insert_remap_to_gud_db(user, host, port, db,
    directory, source_name):

    # Initialize
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

    # Initialize TF-binding table
    if not engine.has_table("tf_binding"):
        raise ValueError("GUD db does not have \"tf_binding\" table!")
    table = TfBinding()
    table.metadata.bind = engine
    table.metadata.create_all(engine)
    mapper(Model, table.__table__)

    # For each file... #
    for file_name in os.listdir(directory):
        # Initialize
        rows = []
        # Skip if not BED file
        if not file_name.endswith(".bed") and not file_name.endswith(".bed.gz"): continue
        # If compressed file...
        if file_name.endswith(".gz"): gz = True
        else: gz = False
        # Get TF name
        m = re.search("^remap\d+\_(\S+)\_all_macs2", file_name)
        tf_name = m.group(1)
        # For each line...
        for line in GUDglobals.parse_tsv_file(os.path.join(directory, file_name), gz=gz):
            # Ignore non-standard chroms, scaffolds, etc.
            m = re.search("^chr(\S+)$", line[0])
            if not m.group(1) in GUDglobals.chroms: continue
            # Get bin
            start = int(line[1])
            end = int(line[2])
            bin = assign_bin(start, end)
            # Get sample
            m = re.search("^.+\..+\.(.+)$", line[3])
            cell_or_tissue = m.group(1)
            # Create model
            model = Model()
            model.bin = bin
            model.chrom = line[0]
            model.start = start
            model.end = end
            model.tf_name = tf_name
            model.cell_or_tissue = cell_or_tissue
            model.experiment_type = "ChIP-seq"
            model.source_name = source_name
            model.date = today
            # Upsert model & commit
            session.merge(model)
            session.commit()

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    # Parse arguments
    args = parse_args()

    # Insert ReMap data to GUD database
    insert_remap_to_gud_db(args.user, args.host, args.port,
        args.db, args.directory, args.source)
