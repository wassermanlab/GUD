#!/usr/bin/env python2.7

import os, sys, re
import argparse
import ConfigParser
from datetime import date
import pybedtools
from sqlalchemy import create_engine
from sqlalchemy.orm import mapper, scoped_session, sessionmaker
from sqlalchemy_utils import database_exists
import warnings

# Append OnTarget module to path
ontarget_path = os.path.join(os.path.dirname(__file__),
    os.pardir, os.pardir, os.pardir)
sys.path.append(ontarget_path)

# Import from OnTarget
from lib import parse_tsv_file
from lib.GUD.ORM.dna_accessibility import DnaAccessibility
from lib.GUD.ORM.enhancer import Enhancer
from lib.GUD.ORM.histone_modification import HistoneModification
from lib.GUD.ORM.tad import Tad
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

    parser.add_argument("file", help="BED file containing the genomic features")

    feats = ["accessibility", "enhancer", "histone", "tad", "tf"]
    parser.add_argument("feat_type", choices=feats, help="Type of genomic feature (i.e. %s)" % ", ".join(feats), metavar="feature")

    parser.add_argument("sample", help="Sample name (e.g. \"lung fibroblasts\", \"GM12878\", etc.)")
    parser.add_argument("exp_type", help="Type of experiment (e.g. \"CAGE\", \"Hi-C\", etc.)")
    parser.add_argument("source", help="Source name (e.g. \"PMID:29449408\")")

    # Optional args
    parser.add_argument("--histone-type", help="Histone type (e.g. \"H3K27ac\", \"H3K4me1\", etc.)")
    parser.add_argument("--restriction-enzyme", help="Restriction enzyme (e.g. \"HindIII\", \"MboI\", etc.)")
    parser.add_argument("--tf-name", help="TF name (e.g. \"FOS\", \"JUN\", \"POLR2A\", etc.)")

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

def insert_bed_to_gud(user, host, port, db,
    file_name, feat_type, cell_or_tissue, experiment_type,
    source_name, histone_type, restriction_enzyme, tf_name):
    """
    """

    # Initialize
    lines = []
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

    # Initialize TF-binding table
    if feat_type == "accessibility":
        if not engine.has_table("dna_accessibility"):
            raise ValueError("GUD db does not have \"dna_accessibility\" table!")
        table = DnaAccessibility()
    if feat_type == "enhancer":
        if not engine.has_table("enhancer"):
            raise ValueError("GUD db does not have \"enhancer\" table!")
        table = Enhancer()
    if feat_type == "histone":
        if not engine.has_table("histone_modification"):
            raise ValueError("GUD db does not have \"histone_modification\" table!")
        if not histone_type:
            raise ValueError("A histone type must be provided!")
        table = HistoneModification()
    if feat_type == "tad":
        if not engine.has_table("tad"):
            raise ValueError("GUD db does not have \"tad\" table!")
        if not restriction_enzyme:
            warnings.warn("\nRestriction enzyme is not provided...\n\tUsing unknown...\n")
            restriction_enzyme = "Unknown"
        hierarchical_level = 0
        table = Tad()
    if feat_type == "tf":
        if not engine.has_table("tf_binding"):
            raise ValueError("GUD db does not have \"tf_binding\" table!")
        if not tf_name:
            raise ValueError("A TF name must be provided!")
        table = TfBinding()
    table.metadata.bind = engine
    table.metadata.create_all(engine)
    mapper(Model, table.__table__)

    # Get lines
    if file_name.endswith(".gz"): gz = True
    else: gz = False
    # For each line...
    for line in parse_tsv_file(file_name, gz=gz):
        lines.append("\t".join(line[:3]))

    # If lines...
    if lines:
        # Create BED object
        bed_obj = pybedtools.BedTool("\n".join(lines), from_string=True)
        # Sort BED object
        for chrom, start, end in bed_obj.sort():
            # Ignore non-standard chroms, scaffolds, etc.
            m = re.search("^chr\w{1,2}$", chrom)
            if not m: continue
            # Get bins
            bins = BinRange().allBinsInRange(int(start), int(end))
            # For each bin...
            for bin in bins:
                # Create model
                model = Model()
                model.bin = bin
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
                    model.hierarchical_level = hierarchical_level
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

    # Insert BED file to GUD
    insert_bed_to_gud(args.user, args.host, args.port,
        args.db, args.file, args.feat_type, args.sample,
        args.exp_type, args.source, args.histone_type,
        args.restriction_enzyme, args.tf_name)
