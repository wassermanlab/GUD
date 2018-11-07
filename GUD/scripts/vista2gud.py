#!/usr/bin/env python

import os, sys, re
import argparse
from binning import assign_bin
from datetime import date
import getpass
import pybedtools
from sqlalchemy import create_engine
from sqlalchemy.orm import mapper, scoped_session, sessionmaker
from sqlalchemy_utils import database_exists
import warnings

# Import from GUD module
from GUD import GUDglobals
from GUD.ORM.enhancer import Enhancer

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

    parser = argparse.ArgumentParser(description="this script inserts \"VISTA\" enhancers into GUD.")

    parser.add_argument("file", help="FASTA file", metavar="file")

    # Optional args
    parser.add_argument("-e", "--exp-type", default="in vivo",
        help="Type of experiment (e.g. \"mouse transgenic enhancer assay\"; default= \"in vivo\")")
    parser.add_argument("--source", default="VISTA",
        help="Source name (e.g. \"PMID:17086198\"; default = \"VISTA\")")

    # MySQL args
    mysql_group = parser.add_argument_group("mysql arguments")
    mysql_group.add_argument("-d", "--db", default="hg19",
        help="Database name (default = \"hg19\")")
    mysql_group.add_argument("-H", "--host", default="localhost",
        help="Host name (default = \"localhost\")")
    mysql_group.add_argument("-P", "--port", default=5506, type=int,
        help="Port number (default = 5506)")
    mysql_group.add_argument("-u", "--user", default=getpass.getuser(),
        help="User name (default = current user)")

    return parser.parse_args()

def insert_vista_to_gud_db(user, host, port, db, fasta_file,
    experiment_type, source_name):

    # Initialize
    lines = []
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

    # Initialize enhancer table
    if not engine.has_table("enhancer"):
        raise ValueError("GUD db does not have \"enhancer\" table!")
    table = Enhancer()
    table.metadata.bind = engine
    table.metadata.create_all(engine)
    mapper(Model, table.__table__)

    # For each header, sequence...
    for header, sequence in GUDglobals.parse_fasta_file(fasta_file):
        # Skip negative enhancers
        m = re.search("|\s+negative\s+|", header)
        if m: print("here")
        else: print("there")
        continue
        # Get chrom, start, end
        m = re.search("Human|chr16:86430087-86430726 | element 1 | positive  | neural tube[12/12] | hindbrain (rhombencephalon)[12/12] | limb[3/12] | cranial nerve[8/12]", header)
        print(header)
        exit(0)
        # Sort BED object
        for chrom, start, end in bed_obj.sort():
            # Create model
            model = Model()
            model.bin = assign_bin(int(start), int(end))
            model.chrom = chrom
            model.start = start
            model.end = end
            model.cell_or_tissue = cell_or_tissue
            model.experiment_type = experiment_type
            model.source_name = source_name
            model.date = today
            if feat_type == "histone":
                model.histone_type = histone_type
            if feat_type == "tad":
                model.restriction_enzyme = restriction_enzyme
            if feat_type == "tf":
                model.tf_name = tf_name
            # Upsert model & commit
            session.merge(model)
            session.commit()

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    # Parse arguments
    args = parse_args()

    # Insert VISTA to GUD database
    insert_vista_to_gud_db(args.user, args.host, args.port,
        args.db, args.file, args.exp_type, args.source)