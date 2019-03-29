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

# Import from GUD
from GUD import GUDglobals
from GUD.ORM.chrom import Chrom
from GUD.ORM.experiment import Experiment
from GUD.ORM.region import Region
from GUD.ORM.sample import Sample
from GUD.ORM.source import Source

usage_msg = """
usage: initialize.py --genome STR [-h] 
                     [-d STR] [-H STR] [-p STR] [-P STR] [-u STR]                            
"""

help_msg = """%s

initializes a GUD database for the given genome.

  --genome STR        genome assembly

optional arguments:
  -h, --help          show this help message and exit

mysql arguments:
  -d STR, --db STR    database name (default = "%s")
  -H STR, --host STR  host name (default = "localhost")
  -p STR, --pwd STR   password (default = ignore this option)
  -P STR, --port STR  port number (default = %s)
  -u STR, --user STR  user name (default = current user)
""" % \
(
    usage_msg,
    GUDglobals.db_name,
    GUDglobals.db_port
)

#-------------#
# Functions   #
#-------------#

def parse_args():
    """
    This function parses arguments provided via
    the command line and returns an {argparse}
    object.
    """

    parser = argparse.ArgumentParser(
        add_help=False,
    )

    # Mandatory arguments
    parser.add_argument("--genome")

    # Optional args
    optional_group = parser.add_argument_group(
        "optional arguments"
    )
    optional_group.add_argument(
        "-h", "--help",
        action="store_true"
    )
    
    # MySQL args
    mysql_group = parser.add_argument_group(
        "mysql arguments"
    )
    mysql_group.add_argument(
        "-d", "--db",
        default=GUDglobals.db_name,
    )
    mysql_group.add_argument(
        "-H", "--host",
        default="localhost"
    )
    mysql_group.add_argument("-p", "--pwd")
    mysql_group.add_argument(
        "-P", "--port",
        default=GUDglobals.db_port
    )
    mysql_group.add_argument(
        "-u", "--user",
        default=getpass.getuser()
    )

    args = parser.parse_args()

    check_args(args)

    return args

def check_args(args):
    """
    This function checks an {argparse} object.
    """

    # Print help
    if args.help:
        print(help_msg)
        exit(0)

    # Check mandatory arguments
    if (not args.genome):
        print(": "\
            .join(
                [
                    "%s\ninitialize.py" % usage_msg,
                    "error",
                    "argument \"--genome\" is required\n"
                ]
            )
        )
        exit(0)

    # Check MySQL password
    if not args.pwd:
        args.pwd = ""

def main():

    # Parse arguments
    args = parse_args()

    # Initialize GUD database
    initialize_gud_db(
        args.user,
        args.pwd,
        args.host,
        args.port,
        args.db,
        args.genome
    )

def initialize_gud_db(user, pwd, host, port,
    db, genome):

    # Initialize
    db_name = "mysql://{}:{}@{}:{}/{}".format(
        user, pwd, host, port, db
    )
    if not database_exists(db_name):
        create_database(db_name)
    session = scoped_session(sessionmaker())
    engine = create_engine(
        db_name,
        echo=False
    )
    session.remove()
    session.configure(
        bind=engine,
        autoflush=False,
        expire_on_commit=False
    )

    table = Chrom()
    if not engine.has_table(
        table.__tablename__
    ):
        # Intialize
        rows = []
        # Create table
        table.__table__.create(bind=engine)
        # Get UCSC FTP file
        directory, file_name =\
            get_ftp_dir_and_file(
                genome,
                "chrom_size"
            )
        # Download data
        for line in fetch_lines_from_ftp_file(
            genome,
            directory,
            file_name
        ):
            # Split line
            line = line.split("\t")
            # Ignore non-standard chroms,
            # scaffolds, etc.
            m = re.search("^chr(\S+)$", line[0])
            if not m.group(1) in GUDglobals.chroms:
                continue
            # Add row
            rows.append(
                {
                    "chrom": line[0],
                    "size": line[1]
                }
            )
        # Insert rows to table
        engine.execute(
            table.__table__.insert(),
            rows
        )

    table = Experiment()
    if not engine.has_table(
        table.__tablename__
    ):
        # Create table
        table.__table__.create(
            bind=engine
        )

    table = Region()
    if not engine.has_table(
        table.__tablename__
    ):
        # Create table
        table.__table__.create(
            bind=engine
        )

    table = Sample()
    if not engine.has_table(
        table.__tablename__
    ):
        # Create table
        table.__table__.create(
            bind=engine
        )

    table = Source()
    if not engine.has_table(
        table.__tablename__
    ):
        # Create table
        table.__table__.create(
            bind=engine
        )

def get_ftp_dir_and_file(genome, data_type):

    # Initialize
    ftp = FTP("hgdownload.soe.ucsc.edu")
    ftp.login()

    # Change into "genome" folder
    try:
        ftp.cwd(
            os.path.join(
                "goldenPath",
                genome
            )
        )
    except:
        print(": "\
            .join(
                [
                    "\ninitialize.py",
                    "error",
                    "cannot connect to UCSC's FTP goldenPath folder",
                    "\"%s\"\n" % genome
                ]
            )
        )
        exit(0)

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
        ftp.cwd(
            os.path.join(
                "goldenPath",
                genome,
                directory
            )
        )
    except:
        print(": "\
            .join(
                [
                    "\ninitialize.py",
                    "error",
                    "cannot connect to UCSC's FTP goldenPath folder",
                    "\"%s/%s\"\n" % (genome, directory)
                ]
            )
        )
        exit(0)

    # If valid file...
    if file_name in ftp.nlst():
        # Retrieve FTP file
        ftp.retrbinary(
            "RETR %s" % file_name,
            callback=handle_bytes
        )
        BIO.seek(0) # Go back to the start
        # If compressed file...
        if file_name.endswith(".gz"):
            f = gzip.GzipFile(
                fileobj=BIO,
                mode="rb"
            )
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

if __name__ == "__main__": main()