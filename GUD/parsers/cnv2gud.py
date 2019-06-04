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
from GUD.ORM.copy_number_variant import CNV
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

    parser.add_argument("cvn_file", help="cnv tsv file")

    # Optional args
    parser.add_argument("--source", default="ISCA-ClinGen", help="Source name")

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


def insert_cnv_to_gud_db(user, host, port, db, tsv_file, source_name):

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
    table = CNV()
    table.metadata.bind = engine
    try:
        table.metadata.create_all(engine)
    except:
        raise ValueError("Cannot create \"copy_number_variants\" table!")
    if not engine.has_table("regions"):
        raise ValueError("No regions table!")
    if not engine.has_table("sources"):
        raise ValueError("No sources table!")

    # get source
    source = Source()
    sou = source.select_by_name(session, source_name)
    if not sou: 
        source.name = source_name
        session.merge(source)
        session.commit()
        sou = source.select_by_name(session, source_name)

    chroms = ["chr"+ch for ch in GUDglobals.chroms]
    # parse table
    with open(tsv_file) as f:
        for line in f:
            if not line.startswith("#"):
                split_line = line.split("\t")
                split_line[-1] = split_line[-1].rstrip()
		        chrom = str(split_line[0]) 
                start = int(split_line[1])
                end = int(split_line[2])
                name = str(split_line[3])
                clinical_interpretation = str(split_line[4])
                variant_type = str(split_line[5])
                copy_number = int(split_line[6])
                 
                if chrom in chroms:
                    # region entry 
                    region = Region()
                    reg = region.select_by_exact_location(session, chrom, start, end)
                    if not reg: 
                        region.bin = assign_bin(start, end)
                        region.chrom = chrom
                        region.start = start
                        region.end = end 
                        session.merge(region)
                        session.commit()
                        reg = region.select_by_exact_location(session, chrom, start, end)
    
                    # str entry 
                    cnv = CNV()
                    if cnv.is_unique(session, name):
                        cnv.uid = name 
                        cnv.copy_number = copy_number
                        cnv.variant_type = variant_type
			            cnv.clinical_interpretation = clinical_interpretation
                        cnv.regionID = reg.uid
                        cnv.sourceID = sou.uid
                        session.merge(cnv)
                        session.commit()

#-------------#
# Main        #
#-------------#


if __name__ == "__main__":

    # Parse arguments
    args = parse_args()

    # Insert ENCODE data to GUD database
    insert_cnv_to_gud_db(args.user, args.host, args.port, args.db,
                         args.cvn_file, args.source)
