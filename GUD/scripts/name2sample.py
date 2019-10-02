#!/usr/bin/env python

import argparse
import os
import shutil

# Import from GUD module
from GUD import GUDUtils
from GUD.ORM.sample import Sample
from GUD.parsers import ParseUtils

usage_msg = """
usage: name2sample.py --name STR [-h] [--dummy-dir DIR] [-o FILE]
                      [-d STR] [-H STR] [-p STR] [-P STR] [-u STR]
"""

help_msg = """%s
identifies one or more samples in GUD that match the given name.

  --name STR          name to match (e.g. "B cell")

optional arguments:
  -h, --help          show this help message and exit
  --dummy-dir DIR     dummy directory (default = "/tmp/")
  -o FILE             output file (default = stdout)

mysql arguments:
  -d STR, --db STR    database name (default = "%s")
  -H STR, --host STR  host name (default = "%s")
  -p STR, --pwd STR   password (default = ignore this option)
  -P STR, --port STR  port number (default = %s)
  -u STR, --user STR  user name (default = "%s")
""" % \
(
    usage_msg,
    GUDUtils.db,
    GUDUtils.host,
    GUDUtils.port,
    GUDUtils.user
)

#-------------#
# Functions   #
#-------------#

def parse_args():
    """
    This function parses arguments provided via the command line and returns
    an {argparse} object.
    """

    parser = argparse.ArgumentParser(add_help=False)

    # Mandatory arguments
    parser.add_argument("--name")

    # Optional args
    optional_group = parser.add_argument_group("optional arguments")
    optional_group.add_argument("-h", "--help", action="store_true")
    optional_group.add_argument("--dummy-dir", default="/tmp/")
    optional_group.add_argument("-o")

    # MySQL args
    mysql_group = parser.add_argument_group("mysql arguments")
    mysql_group.add_argument("-d", "--db", default=GUDUtils.db)
    mysql_group.add_argument("-H", "--host", default=GUDUtils.host)
    mysql_group.add_argument("-p", "--pwd")
    mysql_group.add_argument("-P", "--port", default=GUDUtils.port)
    mysql_group.add_argument("-u", "--user", default=GUDUtils.user)

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
    if not args.name:
        print(": "\
            .join(
                [
                    "%s\nname2sample.py" % usage_msg,
                    "error",
                    "argument \"--name\" is required\n"
                ]
            )
        )
        exit(0)

def main():

    # Parse arguments
    args = parse_args()

    # Initialize
    samples = set()
    dummy_file = os.path.join(
        os.path.abspath(args.dummy_dir),
        "%s.%s.txt" % (os.path.basename(__file__), os.getpid())
    )

    # Set MySQL options
    GUDUtils.user = args.user
    GUDUtils.pwd = args.pwd
    GUDUtils.host = args.host
    GUDUtils.port = args.port
    GUDUtils.db = args.db

    # Get database name
    db_name = GUDUtils._get_db_name()

    # Get engine/session
    engine, Session = GUDUtils._get_engine_session(db_name)

    # Start a new session
    session = Session()

    # For each sample...
    for sample, score in Sample.select_by_fuzzy_string_matching(session, args.name):

        # If sample not in samples...
        if sample not in samples:

            # Write
            ParseUtils.write(dummy_file, "%s\t%s" % (sample.name, score))

        # Add sample
        samples.add(sample)

    if samples:

        # If output file...
        if args.o:
            shutil.copy(dummy_file, os.path.abspath(args.o))
        # ... Else, print on stdout...
        else:
            for line in ParseUtils.parse_file(dummy_file):
                ParseUtils.write(None, line)

        # Delete dummy file
        os.remove(dummy_file)

    else:
        raise ValueError("No samples found!!!")

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()