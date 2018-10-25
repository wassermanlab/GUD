#!/usr/bin/env python2.7

import os, sys, re
import argparse
import ConfigParser
import pybedtools
from sqlalchemy import create_engine
from sqlalchemy.orm import mapper, scoped_session, sessionmaker
from sqlalchemy_utils import database_exists

# Append OnTarget module to path
ontarget_path = os.path.join(os.path.dirname(__file__),
    os.pardir, os.pardir, os.pardir)
sys.path.append(ontarget_path)

# Import from OnTarget
from lib import parse_tsv_file
from lib.GUD.ORM.enhancer import Enhancer
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

    parser.add_argument("file", help="File containing the \"Individual Experiments\" metadata (i.e. https://gro-seq.colorado.edu/data-analysis/exp-downloads-info/)")
    parser.add_argument("dir", help="Directory containing the Tfit calls of \"Individual Experiments\" (e.g. https://nascent.colorado.edu/EMG_out/human/human_EMGs.tgz)")

    # Optional args
    parser.add_argument("--source", default="PMID:29449408", help="Source name (e.g. \"Dowell Lab\"; default = \"PMID:29449408\")")
    parser.add_argument("--date", default="2018-02-15", help="Publication date (in YYYY-MM-DD format; default = \"2018-02-15\")")

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

def insert_groseq_to_gud(user, host, port, db,
    file_name, dir_name, source_name, date_published):
    """
    """

    # Initialize
    models = []
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

    # If enhancer table does not exist...
    if not engine.has_table("enhancer"):
        raise ValueError("GUD db does not have \"enhancer\" table!")

    # For each line...
    for line in parse_tsv_file(file_name):
        srr = line[1]
        cell_or_tissue = line[5]
        if line[4].upper() != cell_or_tissue.upper():
            cell_or_tissue += ", %s" % line[4]
        cell_or_tissue += ", %s" % line[3]
        experiment_type = line[-1]
        metadata.setdefault((cell_or_tissue, experiment_type), [])
        metadata[(cell_or_tissue, experiment_type)].append(srr)

    # Initialize TF-binding table
    table = Enhancer()
    table.metadata.bind = engine
    table.metadata.create_all(engine)
    mapper(Model, table.__table__)
    # For each cell/tissue and experiment...
    for cell_or_tissue, experiment_type in sorted(metadata):
        # Initialize
        lines = []
        # For each SRR...
        for srr in sorted(metadata[(cell_or_tissue, experiment_type)]):                
            # If SRR file exists
            file_name = os.path.join(dir_name, "%s-1_divergent_classifications.bednoinvalidlines.bed" % srr)
            if os.path.exists(file_name):
                # For each line...
                for line in parse_tsv_file(file_name):
                    lines.append("\t".join([line[0], line[1], line[2]]))
        # If lines...
        if lines:
            # Create BED object
            bed_obj = pybedtools.BedTool("\n".join(lines), from_string=True)
            # Sort BED object
            sorted_bed_obj = bed_obj.sort()
            # Merge
            for chrom, start, end in sorted_bed_obj.merge():
                # Ignore non-standard chroms, scaffolds, etc.
                m = re.search("^chr\w{1,2}$", chrom)
                if not m: continue
                # Get bins
                bins = BinRange().allBinsInRange(int(start), int(end))
                # Create model
                model = Model()
                model.bin = bins[0]
                model.chrom = chrom
                model.start = start
                model.end = end
                model.cell_or_tissue = cell_or_tissue
                model.experiment_type = experiment_type
                model.source_name = source_name
                model.date_published = date_published
                # Add model 
                models.append(model)
                # Merge models in bulks of 100,000
                if len(models) == 100000:
                    # For each existing model...
                    for model in models:
                        session.merge(model)
                        # Commit
                        session.commit()
                    # Clear models
                    models = []

    # Merge remaining models
    for model in models:
        session.merge(model)
        # Commit
        session.commit()

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    # Parse arguments
    args = parse_args()

    # Insert GRO-seq data to enhancer table
    insert_groseq_to_gud(args.user, args.host, args.port,
        args.db, args.file, args.dir, args.source, args.date)
