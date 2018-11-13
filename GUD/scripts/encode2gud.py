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
from GUD.ORM.dna_accessibility import DnaAccessibility
from GUD.ORM.histone_modification import HistoneModification
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

    parser = argparse.ArgumentParser(description="this script inserts \"accessibility\", \"histone\" or \"tf\" data from ENCODE into GUD. arguments \"metadata\" and \"directory\" refer to the execution \"xargs -n 1 curl -O -L < file.txt\". genomic features include \"accessibility\", \"histone\" and \"tf\".")

    parser.add_argument("genome", help="Genome assembly")
    parser.add_argument("metadata", help="Metadata file")
    parser.add_argument("directory", help="Downloads directory")

    feats = ["accessibility", "histone", "tf"]
    parser.add_argument("feat_type", choices=feats, help="Type of genomic feature", metavar="feature_type")

    # Optional args
    parser.add_argument("--source", default="ENCODE", help="Source name (e.g. \"PMID:22955616\"; default = \"ENCODE\")")

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

def insert_encode_to_gud_db(user, host, port, db, genome,
    metadata_file, directory, feat_type, source_name):

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
    if feat_type == "accessibility":
        if not engine.has_table("dna_accessibility"):
            raise ValueError("GUD db does not have \"dna_accessibility\" table!")
        table = DnaAccessibility()
    if feat_type == "histone":
        if not engine.has_table("histone_modification"):
            raise ValueError("GUD db does not have \"histone_modification\" table!")
        table = HistoneModification()
    if feat_type == "tf":
        if not engine.has_table("tf_binding"):
            raise ValueError("GUD db does not have \"tf_binding\" table!")
        table = TfBinding()
    table.metadata.bind = engine
    table.metadata.create_all(engine)
    mapper(Model, table.__table__)

    # For each line...
    for line in GUDglobals.parse_tsv_file(metadata_file):
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
        if len(line[9]) > 0 or len(line[10]) > 0 or len(line[11]) > 0:
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
#        # Warn audits
#        if audit:
#            warnings.warn("\nAudition for sample (%s) detected a problem: \"%s\"\n\tConsider skipping...\n" % (
#                accession, audit))
        # This is a released sample!
        if assembly == genome and status == "released":
            metadata.setdefault((cell_or_tissue, experiment_type, experiment_target), [])
            metadata[(cell_or_tissue, experiment_type, experiment_target)].append(accession)

    # For each cell/tissue, experiment and target...
    for cell_or_tissue, experiment_type, experiment_target in sorted(metadata):
        # Initialize
        lines = []
        # For each accession...
        for accession in sorted(metadata[(cell_or_tissue, experiment_type, experiment_target)]):                
            # If accession file exists
            file_name = os.path.join(directory, "%s.bed.gz" % accession)
            if os.path.exists(file_name):
                try:
                    # For each line...
                    for line in GUDglobals.parse_tsv_file(file_name, gz=True):
                        # Skip if not enough elements
                        if len(line) < 3: continue
                        # Ignore non-standard chroms, scaffolds, etc.
                        m = re.search("^chr(\S+)$", line[0])
                        if not m.group(1) in GUDglobals.chroms: continue
                        # Skip if not start or end
                        if not line[1].isdigit(): continue
                        if not line[2].isdigit(): continue
                        # If start is smaller than end
                        if int(line[1]) < int(line[2]):
                            lines.append("\t".join(line[:3]))
                except:
                    warnings.warn("\nCould not read file: \"%s\"\n\tSkipping file...\n" % file_name)
        # If lines...
        if lines:
            # Create BED object
            bed_obj = pybedtools.BedTool("\n".join(lines), from_string=True)
            # Sort and merge
            for chrom, start, end in bed_obj.sort().merge():
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
                    model.histone_type = experiment_target
                if feat_type == "tf":
                    model.tf_name = experiment_target
                # Upsert model & commit
                session.merge(model)
                session.commit()
            # Empty cache
            pybedtools.cleanup()

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    # Parse arguments
    args = parse_args()

    # Insert ENCODE data to GUD database
    insert_encode_to_gud_db(args.user, args.host, args.port,
        args.db, args.genome, args.metadata, args.directory,
        args.feat_type, args.source)