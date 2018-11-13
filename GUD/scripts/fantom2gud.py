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
#from GUD.ORM.tss import TSS

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

    parser = argparse.ArgumentParser(description="this script inserts \"enhancer\" or \"tss\" data from FANATOM into GUD. download \"http://slidebase.binf.ku.dk/human_enhancers/presets/serve/hg19_permissive_enhancers_expression_rle_tpm.csv.gz\" and \"\" for enhancer and tss data, respectively.")

    parser.add_argument("matrix", help="Expression (TPM/RLE normalized) matrix across all FANTOM libraries")

    feats = ["enhancer", "tss"]
    parser.add_argument("feat_type", choices=feats, help="Type of genomic feature", metavar="feature_type")

    # Optional args
    parser.add_argument("-b", "--bed", help="BED file of features on which to focus (e.g. \"robust_enhancers.bed\")")
    parser.add_argument("--source", default="FANTOM", help="Source name (e.g. \"PMID:24670763\" for TSSs or \"PMID:24670764\" for enhancers; default = \"FANTOM\")")

    # MySQL args
    mysql_group = parser.add_argument_group("mysql arguments")
    mysql_group.add_argument("-d", "--db", default="hg19",
        help="Database name (default = \"hg19\")")
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

def insert_fantom_to_gud_db(user, host, port, db, matrix_file,
    feat_type, source_name, bed_file=None):

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
    if matrix_file.endswith(".gz"): gz = True
    else: gz = False

    # Initialize table
    if feat_type == "enhancer":
        if not engine.has_table("enhancer"):
            raise ValueError("GUD db does not have \"enhancer\" table!")
        table = Enhancer()
#    if feat_type == "tss":
#        if not engine.has_table("tss"):
#            raise ValueError("GUD db does not have \"tss\" table!")
#        table = Tss()
    table.metadata.bind = engine
    table.metadata.create_all(engine)
    mapper(Model, table.__table__)

    # If BED file...
    if bed_file:
        # If BED file exists...
        if os.path.exists(bed_file):
            try:
                # Create BED object
                bed_obj = pybedtools.BedTool(bed_file)
            except:
                warnings.warn("\nCould not read file: \"%s\"\n\tSkipping file...\n" % file_name)

    # For each line...
    for line in GUDglobals.parse_csv_file(matrix_file, gz):
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
        print(line)
        exit(0)
#        accession = line[0]
#        experiment_type = line[4]
#        cell_or_tissue = "%s %s" % (line[6], line[7])
#        experiment_target = None
#        m = re.search("^(.+)-(human|mouse)$", line[12])
#        if m: experiment_target = m.group(1)
#        treatment = None
#        if len(line[9]) > 0 or len(line[10]) > 0 or len(line[11]) > 0:
#            treatment = "%s %s %s" % (line[9], line[10], line[11])
#        assembly = line[37]
#        status = line[40]
#        audit = None
#        if len(line[41]) > 0 or len(line[42]) > 0:
#            audit = "%s|%s" % (line[41], line[42])
#        # Skip treated samples
#        if treatment:
#            warnings.warn("\nSample (%s) received treatment: \"%s\"\n\tSkipping sample...\n" % (
#                accession, treatment))
#            continue
##        # Warn audits
##        if audit:
##            warnings.warn("\nAudition for sample (%s) detected a problem: \"%s\"\n\tConsider skipping...\n" % (
##                accession, audit))
#        # This is a released sample!
#        if assembly == genome and status == "released":
#            metadata.setdefault((cell_or_tissue, experiment_type, experiment_target), [])
#            metadata[(cell_or_tissue, experiment_type, experiment_target)].append(accession)
#
#    # For each cell/tissue, experiment and target...
#    for cell_or_tissue, experiment_type, experiment_target in sorted(metadata):
#        # Initialize
#        lines = []
#        # For each accession...
#        for accession in sorted(metadata[(cell_or_tissue, experiment_type, experiment_target)]):                
#
#        # If lines...
#        if lines:
#            # Create BED object
#            bed_obj = pybedtools.BedTool("\n".join(lines), from_string=True)
#            # Sort and merge
#            for chrom, start, end in bed_obj.sort().merge():
#                # Create model
#                model = Model()
#                model.bin = assign_bin(int(start), int(end))
#                model.chrom = chrom
#                model.start = start
#                model.end = end
#                model.cell_or_tissue = cell_or_tissue
#                model.experiment_type = experiment_type
#                model.source_name = source_name
#                model.date = today
#                if feat_type == "histone":
#                    model.histone_type = experiment_target
#                if feat_type == "tf":
#                    model.tf_name = experiment_target
#                # Upsert model & commit
#                session.merge(model)
#                session.commit()
#            # Empty cache
#            pybedtools.cleanup()

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    # Parse arguments
    args = parse_args()

    # Insert FANTOM data to GUD database
    insert_fantom_to_gud_db(args.user, args.host, args.port,
        args.db, args.matrix, args.feat_type, args.source, args.bed)