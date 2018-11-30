#!/usr/bin/env python
import os, sys, re
import argparse
from binning import assign_bin
import getpass
import pybedtools
from sqlalchemy import create_engine
from sqlalchemy.orm import mapper, scoped_session, sessionmaker
from sqlalchemy_utils import database_exists
import warnings

# Import from GUD module
from GUD import GUDglobals
from GUD2.ORM.region import Region
from GUD2.ORM.sample import Sample
from GUD2.ORM.source import Source
from GUD2.ORM.experiment import Experiment
from GUD2.ORM.dna_accessibility import DNAAccessibility
from GUD2.ORM.enhancer import Enhancer
from GUD2.ORM.histone_modification import HistoneModification
from GUD2.ORM.tad import TAD
from GUD2.ORM.tf_binding import TFBinding

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

    parser = argparse.ArgumentParser(description="this script inserts \"accessibility\", \"enhancer\", \"histone\", \"tad\" or \"tf\" data from input BED file into GUD. genomic features include \"accessibility\", \"enhancer\", \"histone\", \"tad\" and \"tf\".")

    parser.add_argument("files", nargs="*", help="BED file(s)", metavar="file")

    feats = ["accessibility", "enhancer", "histone", "tad", "tf"]
    parser.add_argument("feat_type", choices=feats, help="Type of genomic feature",
        metavar="feature_type")

    parser.add_argument("sample", help="Sample name (e.g. \"lung fibroblasts\")")
    parser.add_argument("exp_type", help="Type of experiment (e.g. \"CAGE\")")
    parser.add_argument("source", help="Source name (e.g. \"PMID:29449408\")")

    # Optional args
    parser.add_argument("--histone", help="Histone type (e.g. \"H3K27ac\")")
    parser.add_argument("--enzyme", help="Restriction enzyme (e.g. \"HindIII\")")
    parser.add_argument("--tf-name", help="TF name (e.g. \"FOS\")")
    parser.add_argument("--merge", action="store_true", 
        help="Merge overlapping or \"book-ended\" features (default = False)")

    # MySQL args
    mysql_group = parser.add_argument_group("mysql arguments")
    mysql_group.add_argument("-d", "--db", default="hg19",
        help="Database name (default = \"hg19\")")
    mysql_group.add_argument("-H", "--host", default="localhost",
        help="Host name (default = \"localhost\")")
    mysql_group.add_argument("-P", "--port", default=5506, type=int,
        help="Port number (default = 5506)")
    mysql_group.add_argument("-u", "--user", default=getpass.getuser(),
        help="User name (default = current user)")

    args = parser.parse_args()

    if args.feat_type == "histone" and not args.histone:
        raise ValueError("A histone type must be provided!")

    if args.feat_type == "tad" and not args.enzyme:
        warnings.warn("\nA restriction enzyme was not provided...\n")
        warnings.warn("\nSetting \"restriction_enzyme\" field to \"Unknown\"...\n")
        args.enzyme = "Unknown"

    if args.feat_type == "tf" and not args.tf_name:
        raise ValueError("A TF name must be provided!")

    return args

def insert_bed_to_gud_db(user, host, port, db, bed_files,
    feat_type, cell_or_tissue, experiment_type, source_name,
    histone_type=None, restriction_enzyme=None, tf_name=None,
    merge=False):

    # Initialize
    lines = []
    db_name = "mysql://{}:@{}:{}/{}".format(
        user, host, port, db)
    if not database_exists(db_name):
        raise ValueError("GUD db does not exist: %s" % db_name)
    session = scoped_session(sessionmaker())
    engine = create_engine(db_name, echo=False)
    session.remove()
    session.configure(bind=engine, autoflush=False,
        expire_on_commit=False)

    # Initialize table
    if feat_type == "accessibility":
        table = DNAAccessibility()
    if feat_type == "enhancer":
        table = Enhancer()
    if feat_type == "histone":
        if not histone_type:
            raise ValueError("A histone type must be provided!")
        table = HistoneModification()
    if feat_type == "tad":
        if not restriction_enzyme:
            warnings.warn("\nA restriction enzyme was not provided...\n")
            warnings.warn("\nSetting \"restriction_enzyme\" field to \"Unknown\"...\n")
            restriction_enzyme = "Unknown"
        table = TAD()
    if feat_type == "tf":
        if not tf_name:
            raise ValueError("A TF name must be provided!")
        table = TFBinding()
    if not engine.has_table(table.__tablename__):
        try:
            table.metadata.bind = engine
            table.metadata.create_all(engine)
        except:
            raise ValueError("Cannot create table: %s" % table.__tablename__)
    mapper(Model, table.__table__)

    # For each BED file...
    for bed_file in bed_files:
        # Get lines
        if bed_file.endswith(".gz"): gz = True
        else: gz = False
        # For each line...
        for line in GUDglobals.parse_tsv_file(bed_file, gz=gz):
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
    # If lines...
    if lines:
        # Create BED object
        bed_obj = pybedtools.BedTool("\n".join(lines), from_string=True)
        # If merge...
        if merge:
            bed_obj = bed_obj.sort().merge()
        # Sort BED object
        for chrom, start, end in bed_obj.sort():
            # Region
            region = Region()
            reg = region.select_by_pos(session, chrom, start, end)
            if not reg: 
                region.bin = assign_bin(start, end)
                region.chrom = chrom
                region.start = start
                region.end = end 
                session.add(region) 
            # Source
            source = Source()
            sou = source.select_by_name(session, source_name)
            if not sou: 
                source.name = source_name
                session.add(source)
            # Sample
            sample = Sample()
            samp = sample.select_exact_sample(session, name, treatment, cell_line, cancer)
            if not samp: 
                sample.name = cell_or_tissue
                sample.treatment = 
                sample.cell_line = 
                sample.cancer = 
                session.add(sample)
            # Experiment 
            experiment = Experiment()
            exp = experiment.select_by_name(session, name)
            if not exp: 
                experiment.name = experiment_type
                session.add(experiment)
            ## commit region, source, sample, experiment 
            session.commit()
            ## uids
            reg = region.select_by_pos(session, chrom, start, end).uid
            sou = source.select_by_name(session, source_name).uid
            samp = sample.select_exact_sample(session, name, treatment, cell_line, cancer).uid
            exp = experiment.select_by_name(session, name).uid
            ## feature creation 
            model = Model()
            model.regionID = reg
            model.sourceID = sou
            model.sampleID = samp
            model.experimentID = exp
            if feat_type == "accessibility":
                exists = DNAAccessibility().select_unique(session, reg, sou, samp, exp)
                if not exists: 
                    session.add(model)
                    session.commit()
            if feat_type == "enhancer":
                exists = Enhancer().select_unique(session, reg, sou, samp, exp)
                if not exists: 
                    session.add(model)
                    session.commit()
            if feat_type == "histone":
                exists = HistoneModification().select_unique(session, reg, sou, samp, exp, histone_type)
                if not exists: 
                    model.histone_type = histone_type
                    session.add(model)
                    session.commit()
            if feat_type == "tad":
                exists = TAD().select_unique(session, reg, sou, samp, exp, restriction_enzyme)
                if not exists: 
                    model.restriction_enzyme = restriction_enzyme
                    session.add(model)
                    session.commit()
            if feat_type == "tf": 
                exists = TFBinding().select_unique(session, reg, sou, samp, exp, tf_name)
                model.tf = tf_name
                if not exists: 
                    session.add(model)
                    session.commit()

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    # Parse arguments
    args = parse_args()

    # Insert BED file to GUD database
    insert_bed_to_gud_db(args.user, args.host, args.port,
        args.db, args.files, args.feat_type, args.sample,
        args.exp_type, args.source, args.histone,
        args.enzyme, args.tf_name, args.merge)