#!/usr/bin/env python

import os, sys, re
import argparse
from binning import assign_bin
from datetime import date
import getpass
import pybedtools
import shutil
from sqlalchemy import create_engine
from sqlalchemy.orm import mapper, scoped_session, sessionmaker
from sqlalchemy_utils import create_database, database_exists
import warnings

# Import from GUD module
from GUD2 import GUDglobals
from GUD2.ORM.dna_accessibility import DNAAccessibility
from GUD2.ORM.histone_modification import HistoneModification
from GUD2.ORM.region import Region
from GUD2.ORM.sample import Sample
from GUD2.ORM.source import Source
from GUD2.ORM.tf_binding import TFBinding

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
    parser.add_argument("samples", help="ENCODE samples (manually-curated)")

    feats = ["accessibility", "histone", "tf"]
    parser.add_argument("feat_type", choices=feats, metavar="feature",
        help="Type of genomic feature")

    # Optional args
    parser.add_argument("--dummy-dir", default="/tmp/",
        help="Dummy directory (default = /tmp/)")
    parser.add_argument("--source", default="ENCODE",
        help="Source name (e.g. \"PMID:22955616\"; default = \"ENCODE\")")

    # MySQL args
    mysql_group = parser.add_argument_group("mysql arguments")
    mysql_group.add_argument("-d", "--db",
        help="Database name (default = given genome assembly)")
    mysql_group.add_argument("-H", "--host", default="localhost",
        help="Host name (default = localhost)")
    mysql_group.add_argument("-P", "--port", default=5506, type=int,
        help="Port number (default = 5506)")

    user = getpass.getuser()
    mysql_group.add_argument("-u", "--user", default=user,
        help="User name (default = current user)")

    args = parser.parse_args()

    # Set default
    if not args.db:
        args.db = args.genome

    return args

def main():

    # Parse arguments
    args = parse_args()

    # Insert ENCODE data to GUD database
    insert_encode_to_gud_db(args.user, args.host, args.port,
        args.db, args.genome, args.metadata, args.directory,
        args.samples, args.feat_type, args.dummy_dir, args.source)

def insert_encode_to_gud_db(user, host, port, db, genome,
    metadata_file, directory, samples_file, feat_type, dummy_dir,
    source_name):

    # Initialize
    samples = {}
    metadata = {}
    db_name = "mysql://{}:@{}:{}/{}".format(
        user, host, port, db)
    if not database_exists(db_name):
        create_database(db_name)
    session = scoped_session(sessionmaker())
    engine = create_engine(db_name, echo=False)
    session.remove()
    session.configure(bind=engine, autoflush=False,
        expire_on_commit=False)
    today = str(date.today())
    dummy_dir = os.path.join(dummy_dir,
        "%s.%s" % (os.path.basename(__file__), os.getpid()))

    # Get source
    source = Source()
    sou = source.select_by_name(session, source_name)
    if not sou:    
        source.name = source_name
        session.add(source)
        session.commit()
        sou = source.select_by_name(session, source_name)

    # Get samples
    for line in GUDglobals.parse_tsv_file(samples_file):
        # Initialize 
        sample_name = line[0]
        sample_type = line[1]
        cell_or_tissue = line[2]
        if line[3] == "Yes": add = True
        else: add = False
        if line[4] == "Yes": treatment = True
        else: treatment = False
        if line[5] == "Yes": cell_line = True
        else: cell_line = False
        if line[6] == "Yes": cancer = True
        else: cancer = False
        # Get sample
        samples.setdefault(sample_name, {
            "sample_type": sample_type,
            "cell_or_tissue": cell_or_tissue,
            "add": add,
            "treatment": treatment,
            "cell_line": cell_line,
            "cancer": cancer
        })

#    # Initialize table
#    if feat_type == "accessibility":
#        table = DnaAccessibility()
#    if feat_type == "histone":
#        table = HistoneModification()
#    if feat_type == "tf":
#        table = TfBinding()
#    if not engine.has_table(table.__tablename__):
#        try:
#            table.metadata.bind = engine
#            table.metadata.create_all(engine)
#        except:
#            raise ValueError("Cannot create table: %s" % table.__tablename__)
#    mapper(Model, table.__table__)

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
        biosample = line[6]
        experiment_target = None
        m = re.search("^(.+)-(human|mouse)$", line[12])
        if m: experiment_target = m.group(1)
        treatment = None
        if len(line[9]) > 0 or len(line[10]) > 0 or len(line[11]) > 0:
            treatment = "%s %s %s" % (line[9], line[10], line[11])
        assembly = line[37]
        status = line[40]
#        audit = None
#        if len(line[41]) > 0 or len(line[42]) > 0:
#            audit = "%s|%s" % (line[41], line[42])
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
            # Get metadata
            metadata.setdefault((experiment_type, experiment_target), [])
            metadata[(experiment_type, experiment_target)].append((accession, biosample))

    # For each cell/tissue, experiment and target...
    for experiment_type, experiment_target in sorted(metadata):
        if os.path.isdir(dummy_dir): shutil.rmtree(dummy_dir)
        os.mkdir(dummy_dir)
        if experiment_type == "FAIRE-seq":
            for accession, biosample in metadata[(experiment_type, experiment_target)]:
                # Copy BED file
                bed_obj = pybedtools.BedTool(os.path.join(directory, "%s.bed.gz" % accession))
                bed_obj.saveas(os.path.join(dummy_dir, "%s.bed" % accession), compressed=False)
            # Empty cache
            pybedtools.cleanup()
    exit(0)

    # For each cell/tissue, experiment and target...
    for experiment_type, experiment_target in sorted(metadata):
        # Initialize
        lines = []
        merged_lines = []
        # For each accession...
        for accession, biosample in sorted(metadata[(experiment_type, experiment_target)]):                
            # If accession file exists
            file_name = os.path.join(directory, "%s.bed.gz" % accession)
            if os.path.exists(file_name):
                try:
                    # For each line...
                    for line in GUDglobals.parse_tsv_file(file_name, gz=True):
                        # Skip if not enough elements
                        if len(line) < 3: continue
                        # Ignore non-standard chroms, scaffolds, etc.
                        m = re.search("^chr(\w{1,2})$", line[0])
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
                merged_lines.append((chrom, start, end))
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
            # Empty cache
            pybedtools.cleanup()
        print(len(lines), len(merged_lines))

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()