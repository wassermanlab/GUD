#!/usr/bin/env python2.7

import os, sys, re
import argparse
import ConfigParser
import pybedtools
from sqlalchemy import create_engine
from sqlalchemy.orm import scoped_session, sessionmaker
from sqlalchemy_utils import database_exists
import warnings

# Append OnTarget module to path
ontarget_path = os.path.join(os.path.dirname(__file__),
    os.pardir, os.pardir, os.pardir)
sys.path.append(ontarget_path)

# Import from OnTarget
from lib import parse_tsv_file
from lib.GUD.ORM.dna_accessibility import DnaAccessibility
from lib.GUD.ORM.histone_modification import HistoneModification
from lib.GUD.ORM.tf_binding import TfBinding
from lib.GUD.utils.bin_range import BinRange

# Read configuration file
config = ConfigParser.ConfigParser()
config_file = os.path.join(ontarget_path, "config.ini")
config.read(config_file)

#-------------#
# Functions   #
#-------------#

def parse_args():
    """
    This function parses arguments provided via command
    line and returns an {argparse} object.
    """

    parser = argparse.ArgumentParser(description="describe what the script does...")

    parser.add_argument("genome", help="Genome assembly (e.g. \"mm10\")")
    parser.add_argument("file", help="File containing the metadata (i.e. from \"xargs -n 1 curl -O -L < histone-modification.txt\")")
    parser.add_argument("dir", help="Directory containing downloaded files")

    feats = ["accessibility", "histone", "tf"]
    parser.add_argument("feat_type", choices=feats, help="Type of genomic feature (i.e. %s)" % ", ".join(feats), metavar="feature")

    # Optional args
    parser.add_argument("--source", default="ENCODE", help="Source name (e.g. \"PMID:22955616\"; default = \"ENCODE\")")
    parser.add_argument("--date", default="2012-09-06", help="Publication date (in YYYY-MM-DD format; default = \"2012-09-06\")")

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

def insert_encode_to_gud(user, host, port, db, genome,
    file_name, dir_name, feat_type, source_name, date_published):
    """
    """

    # Initialize
    rows = []
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

    # If feature type table does not exist...
    if feat_type == "accessibility":
        if not engine.has_table("dna_accessibility"):
            raise ValueError("GUD db does not have \"dna_accessibility\" table!")
    if feat_type == "histone":
        if not engine.has_table("histone_modification"):
            raise ValueError("GUD db does not have \"histone_modification\" table!")
    if feat_type == "tf":
        if not engine.has_table("tf_binding"):
            raise ValueError("GUD db does not have \"tf_binding\" table!")

    # For each line...
    for line in parse_tsv_file(file_name):
#        line = ['File accession', 'File format', 'Output type', 'Experiment accession', 'Assay',
#                'Biosample term id', 'Biosample term name', 'Biosample type', 'Biosample organism',
#                'Biosample treatments', 'Biosample treatments amount', 'Biosample treatments duration',
#                'Experiment target', 'Library made from', 'Library depleted in', 'Library extraction method',
#                'Library lysis method', 'Library crosslinking method', 'Library strand specific',
#                'Experiment date released', 'Project', 'RBNS protein concentration', 'Library fragmentation method',
#                'Library size range', 'Biological replicate(s)', 'Technical replicate', 'Read length',
#                'Mapped read length', 'Run type', 'Paired end', 'Paired with', 'Derived from', 'Size',
#                'Lab', 'md5sum', 'dbxrefs', 'File download URL', 'Assembly', 'Platform', 'Controlled by',
#                'File Status', 'Audit WARNING', 'Audit INTERNAL_ACTION', 'Audit NOT_COMPLIANT', 'Audit ERROR']
        accession = line[0]
        experiment_type = line[4]
        cell_or_tissue = "%s %s" % (line[6], line[7])
        experiment_target = None
        m = re.search("^(.+)-(human|mouse)$", line[12])
        if m: experiment_target = m.group(1)
        treatment = None
        if len(line[9]) > 0 or len(line[9]) > 0 or len(line[9]) > 0:
            treatment = "%s %s %s" % (line[9], line[10], line[11])
        assembly = line[37]
        status = line[40]
        audit = None
        if len(line[41]) > 0 or len(line[42]) > 0:
            audit = "%s|%s" % (line[41], line[42])
        # Skip treated samples
        if treatment:
            warnings.warn("\nSample (%s) received treatment: \"%s\"\n\tSkipping sample...\n" % (
                accession, treatment))
            continue
        # Warn audits
        if audit:
            warnings.warn("\nAudition for sample (%s) detected a problem: \"%s\"\n\tSkipping sample...\n" % (
                accession, audit))
        # This is a released sample!
        if assembly == genome and status == "released":
            metadata.setdefault((cell_or_tissue, experiment_type, experiment_target), [])
            metadata[(cell_or_tissue, experiment_type, experiment_target)].append(accession)

    # Initialize TF-binding table
    if feat_type == "accessibility": table = DnaAccessibility()
    if feat_type == "histone": table = HistoneModification()
    if feat_type == "tf": table = TfBinding()
    table.metadata.bind = engine
    table.metadata.create_all(engine)
    # For each cell/tissue, experiment and target...
    for cell_or_tissue, experiment_type, experiment_target in sorted(metadata):
        # Initialize
        lines = []
        # For each accession...
        for accession in sorted(metadata[(cell_or_tissue, experiment_type, experiment_target)]):                
            # If accession file exists
            file_name = os.path.join(dir_name, "%s.bed.gz" % accession)
            if os.path.exists(file_name):
                try:
                    # For each line...
                    for line in parse_tsv_file(file_name, gz=True):
                        if int(line[1]) < int(line[2]):
                            lines.append("\t".join([line[0], line[1], line[2]]))
                except:
                    warnings.warn("\nCould not read file: \"%s\"\n\tSkipping file...\n" % file_name)
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
                # Initialize row
                row = {
                    "bin": bins[0],
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "cell_or_tissue": cell_or_tissue,
                    "experiment_type": experiment_type,
                    "source_name": source_name,
                    "date_published": date_published 
                }
                if feat_type == "accessibility":
                    feat_exists = table.feature_exists(
                        session, chrom, start, end, cell_or_tissue,
                        experiment_type, source_name
                    )
                if feat_type == "histone":
                    row.setdefault("histone_type", experiment_target)
                    feat_exists = table.feature_exists(
                        session, chrom, start, end, experiment_target,
                        cell_or_tissue, experiment_type, source_name
                    )
                if feat_type == "tf":
                    row.setdefault("tf_name", experiment_target)
                    feat_exists = table.feature_exists(
                        session, chrom, start, end, experiment_target,
                        cell_or_tissue, experiment_type, source_name
                    )
                # Add row
                if not feat_exists:
                    rows.append(row)
                # Insert rows in bulks of 100,000
                if len(rows) == 100000:
                    engine.execute(table.__table__.insert(), rows)
                    # Clear rows
                    rows = []

    # Insert remaining rows
    engine.execute(table.__table__.insert(), rows)

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    # Parse arguments
    args = parse_args()

    # Insert ENCODE data to GUD
    insert_encode_to_gud(args.user, args.host, args.port,
        args.db, args.genome, args.file, args.dir,
        args.feat_type, args.source, args.date)
