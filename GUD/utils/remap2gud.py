#!/usr/bin/env python2.7

import os, sys, re
import argparse
import ConfigParser
from datetime import date
from sqlalchemy import create_engine
from sqlalchemy.orm import mapper, scoped_session, sessionmaker
from sqlalchemy_utils import database_exists

# Append OnTarget module to path
ontarget_path = os.path.join(os.path.dirname(__file__),
    os.pardir, os.pardir, os.pardir)
sys.path.append(ontarget_path)

# Import from OnTarget
from lib import parse_tsv_file
from lib.GUD.ORM.tf_binding import TfBinding
from lib.GUD.utils.bin_range import BinRange

# Read configuration file
config = ConfigParser.ConfigParser()
config_file = os.path.join(ontarget_path, "config.ini")
config.read(config_file)

#-------------#
# Classes     #
#-------------#

class Model(object): pass

#-------------#
# Functions   #
#-------------#

def parse_args():
    """
    This function parses arguments provided via command
    line and returns an {argparse} object.
    """

    parser = argparse.ArgumentParser(description="describe what the script does...")

    parser.add_argument("dir", help="Directory containing downloaded files (i.e. one BED file per TF)")

    # Optional args
    parser.add_argument("--source", default="ReMap", help="Source name (e.g. \"PMID:29126285\"; default = \"ReMap\")")

    # MySQL args
    mysql_group = parser.add_argument_group("mysql arguments")
    mysql_group.add_argument("-d", "--db", default=config.get("MySQL", "db"),
        help="Database name (e.g. \"mm10\"; default from \"config.ini\" = %s)" % config.get("MySQL", "db"))
    mysql_group.add_argument("-H", "--host", default=config.get("MySQL", "host"),
        help="Host name (e.g. \"ontarget.cmmt.ubc.ca\"; default from \"config.ini\" = %s)" % config.get("MySQL", "host"))
#    mysql_group.add_argument("-p", "--pass", default="", help="User pass")
    mysql_group.add_argument("-P", "--port", default=config.get("MySQL", "port"),
        help="User name (e.g. \"5506\"; default from \"config.ini\" = %s)" % config.get("MySQL", "port"))
    mysql_group.add_argument("-u", "--user", default=config.get("MySQL", "user"),
        help="User name (e.g. \"ontarget_r\"; default from \"config.ini\" = %s)" % config.get("MySQL", "user"))

    return parser.parse_args()

def insert_remap_to_gud(user, host, port, db,
    directory, source_name):
    """
    """

    # Initialize
    merged_models = []
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
        for line in parse_tsv_file(
            os.path.join(directory, file_name), gz=gz):
            # Ignore non-standard chroms, scaffolds, etc.
            m = re.search("^chr\w{1,2}$", line[0])
            if not m: continue
            # Get bins
            start = int(line[1])
            end = int(line[2])
            bins = BinRange().allBinsInRange(start, end)
            # Get sample
            m = re.search("^.+\..+\.(.+)$", line[3])
            cell_or_tissue = m.group(1)
            # For each bin...
            for bin in bins:
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

    # Insert ReMap data to TF-binding table
    insert_remap_to_gud(args.user, args.host, args.port,
        args.db, args.dir, args.source)
