#!/usr/bin/env python

import argparse
from binning import assign_bin
import copy
import getpass
import pybedtools
import os
import re
import shutil
from sqlalchemy import create_engine
from sqlalchemy.orm import (
    scoped_session,
    sessionmaker
)
from sqlalchemy_utils import database_exists
import subprocess
import warnings

# Import from GUD module
from GUD import GUDglobals
from GUD.ORM.dna_accessibility import DNAAccessibility
from GUD.ORM.experiment import Experiment
from GUD.ORM.histone_modification import HistoneModification
from GUD.ORM.region import Region
from GUD.ORM.sample import Sample
from GUD.ORM.source import Source
#from GUD.ORM.tf_binding import TFBinding

#-------------#
# Functions   #
#-------------#

def parse_args():
    """
    This function parses arguments provided via the command
    line using argparse.
    """

    parser = argparse.ArgumentParser(description="inserts \"accessibility\", \"histone\" or \"tf\" data from the ENCODE consortium into GUD. arguments \"metadata\" and \"directory\" refer to the execution \"xargs -n 1 curl -O -L < file.txt\". genomic features include \"accessibility\", \"histone\" and \"tf\".")

    parser.add_argument("genome", help="genome assembly")
    parser.add_argument("metadata", help="metadata file")
    parser.add_argument("directory", help="downloads directory")
    parser.add_argument("samples", help="ENCODE samples (manually-curated)")

    feats = ["accessibility", "histone", "tf"]
    parser.add_argument("feat_type", choices=feats,
        help="type of genomic feature", metavar="feature_type")

    # Optional args
    parser.add_argument("-c", "--cluster", action="store_true",
        help="cluster regions with regCluster (default = False)")
    parser.add_argument("--dummy-dir", default="/tmp/",
        metavar="DUMMYDIR", help="dummy directory (default = /tmp/)")
    parser.add_argument("--source", default="ENCODE",
        help="source name (e.g. \"PMID:22955616\"; default = \"ENCODE\")")

    # MySQL args
    mysql_group = parser.add_argument_group("mysql arguments")
    mysql_group.add_argument("-d", "--db", default="hg19",
        help="database name (default = \"hg19\")")
    mysql_group.add_argument("-H", "--host", default="localhost",
        help="host name (default = localhost)")
    mysql_group.add_argument("-p", "--passwd",
        help="password (default = ignore this option)")
    mysql_group.add_argument("-P", "--port", default=5506, type=int,
        help="port number (default = 5506)")

    user = getpass.getuser()
    mysql_group.add_argument("-u", "--user", default=user,
        help="user name (default = current user)")

    args = parser.parse_args()

    # Set default
    if not args.db:
        args.db = args.genome
    if not args.passwd:
        args.passwd = ""

    return args

def main():

    # Parse arguments
    args = parse_args()

    # Insert ENCODE data to GUD database
    insert_encode_to_gud_db(args.user, args.passwd,
        args.host, args.port, args.db, args.genome,
        args.metadata, args.directory, args.samples,
        args.feat_type, args.cluster, args.dummy_dir,
        args.source)

def insert_encode_to_gud_db(user, passwd, host, port, db,
    genome, metadata_file, directory, samples_file,
    feat_type, cluster, dummy_dir, source_name):

    # Initialize
    samples = {}
    metadata = {}
    db_name = "mysql://{}:{}@{}:{}/{}".format(
        user, passwd, host, port, db)
    if not database_exists(db_name):
        raise ValueError("GUD database does not exist!!!\n\t%s" % db_name)
    session = scoped_session(sessionmaker())
    engine = create_engine(db_name, echo=False)
    session.remove()
    session.configure(bind=engine, autoflush=False,
        expire_on_commit=False)

    # Get source
    source = Source()
    if source.is_unique(session, source_name):    
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

    # Create table
    if feat_type == "accessibility":
        table = DNAAccessibility()
    if feat_type == "histone":
        table = HistoneModification()
    if feat_type == "tf":
        table = TFBinding()
    if not engine.has_table(table.__tablename__):
        # Create table
        table.__table__.create(bind=engine)

    # For each line...
    for line in GUDglobals.parse_tsv_file(metadata_file):
        # Initialize
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
        # Skip treated samples
        if treatment:
            warnings.warn("Skipping treated sample %s (%s)\n" % \
                (accession, treatment)
            )
            continue
#        # Warn audits
#        if audit:
#            warnings.warn("\nAudition for sample (%s) detected a problem: \"%s\"\n\tConsider skipping...\n" % (
#                accession, audit))
        # This is a released sample!
        if assembly == genome and status == "released":
            # Skip sample
            if not samples[biosample]["add"]: continue
            # Get metadata
            if os.path.exists(os.path.join(directory, "%s.bed.gz" % accession)):
                metadata.setdefault((experiment_type, experiment_target), [])
                metadata[(experiment_type, experiment_target)].append((accession, biosample))

    # For each cell/tissue, experiment, target...
    for experiment_type, experiment_target in sorted(metadata):
        # Initialize
        exp_dummy_dir = os.path.join(dummy_dir,
            "%s.%s" % (experiment_type.replace(" ", "_"), experiment_target))
#        # Remove dummy dir
#        if os.path.isdir(exp_dummy_dir): shutil.rmtree(exp_dummy_dir)
        # Create dummy dir
        if not os.path.isdir(exp_dummy_dir): os.mkdir(exp_dummy_dir)
        # Get experiment
        experiment = Experiment()
        if experiment.is_unique(session, experiment_type):
            experiment.name = experiment_type
            session.add(experiment)
            session.commit()
        exp = experiment.select_by_name(session, experiment_type)
        # For each accession, biosample...
        for accession, biosample in metadata[(experiment_type, experiment_target)]:
            # Copy BED file
            gz_bed_file = os.path.join(directory, "%s.bed.gz" % accession)
            bed_file = os.path.join(exp_dummy_dir, "%s.bed" % accession)
            if not os.path.exists(bed_file):
                os.system("zcat %s | sort -k 1,1 -k2,2n > %s" % (gz_bed_file, bed_file))
        # Cluster regions
        if cluster:
            # Initialize
            accession2sample = {}
            label2accession = {}
            regions = []
            bed_files = os.path.join(exp_dummy_dir, "files.txt")
            table_file = os.path.join(exp_dummy_dir, "tableOfTables.txt")
            cluster_file = os.path.join(exp_dummy_dir, "regCluster")
            # Create BED file list
            if not os.path.exists(bed_files):
                # For each file...
                for bed_file in os.listdir(exp_dummy_dir):
                    # Skip non-BED files
                    if not bed_file.endswith(".bed"): continue
                    # Add file to list
                    GUDglobals.write(bed_files, os.path.join(exp_dummy_dir, bed_file))
            # Make table of tables
            if not os.path.exists(table_file):
                process = subprocess.check_output(["regClusterMakeTableOfTables",
                    "uw01", bed_files, table_file], stderr=subprocess.STDOUT)
            # Make clusters
            if not os.path.exists("%s.cluster" % cluster_file):
                process = subprocess.check_output(["regCluster", table_file,
                    "%s.cluster" % cluster_file, "%s.bed" % cluster_file],
                    stderr=subprocess.STDOUT)
            # For each accession, biosample...
            for accession, biosample in metadata[(experiment_type, experiment_target)]:
                # Get sample
                sample = Sample()
                if sample.is_unique(
                    session,
                    samples[biosample]["cell_or_tissue"],
                    samples[biosample]["treatment"],
                    samples[biosample]["cell_line"],
                    samples[biosample]["cancer"]
                ):
                    sample.name = samples[biosample]["cell_or_tissue"]
                    sample.treatment = samples[biosample]["treatment"]
                    sample.cell_line = samples[biosample]["cell_line"]
                    sample.cancer = samples[biosample]["cancer"]
                    session.add(sample)
                    session.commit()
                sam = sample.select_unique(session,
                    samples[biosample]["cell_or_tissue"],
                    samples[biosample]["treatment"],
                    samples[biosample]["cell_line"],
                    samples[biosample]["cancer"])
                accession2sample.setdefault(accession, sam.uid)
            # For each line...
            for line in GUDglobals.parse_tsv_file(table_file):
                m = re.search("%s.%s/(\w+).bed" % (
                    experiment_type, experiment_target), line[0])
                if m: label2accession.setdefault(line[-1], m.group(1))
            # Load BED file
            bed_obj = pybedtools.BedTool("%s.bed" % cluster_file)
            # For each feature...
            for feature in bed_obj:
                # Ignore non-standard chroms, scaffolds, etc.
                m = re.search("^chr(\S+)$", feature[0])
                if not m.group(1) in GUDglobals.chroms: continue
                # Get coordinates
                chrom = feature[0]
                start = int(feature[1])
                end = int(feature[2])
                # Ignore non-standard chroms, scaffolds, etc.
                m = re.search("^chr(\S+)$", chrom)
                if not m.group(1) in GUDglobals.chroms: continue
                # Get region
                region = Region()
                if region.is_unique(session, chrom, start, end):
                    # Insert region
                    region.bin = assign_bin(start, end)
                    region.chrom = chrom
                    region.start = start
                    region.end = end
                    session.add(region)
                    session.commit()
                reg = region.select_unique(session, chrom, start, end)
                regions.append(reg.uid)
            # For each line...
            for line in GUDglobals.parse_tsv_file("%s.cluster" % cluster_file):
                # Get region
                reg_uid = regions[int(line[0]) - 1] 
                # Get sample
                sam_uid = accession2sample[label2accession[line[-1]]]
                # Insert feature
                if feat_type == "accessibility":
                    feat = DNAAccessibility()
                    is_unique = feat.is_unique(
                        session,
                        reg_uid,
                        sam_uid,
                        exp.uid,
                        sou.uid
                    )
#                if feat_type == "histone":
#                    feat = HistoneModification()
#                    is_unique = feat.is_unique(
#                        session,
#                        reg_uid,
#                        sam_uid,
#                        exp.uid,
#                        sou.uid,
#                        experiment_target
#                    )
#                if feat_type == "tf":
#                    feat = TFBinding()
#                    is_unique = feat.is_unique(
#                        session,
#                        reg_uid,
#                        sam_uid,
#                        exp.uid,
#                        sou.uid,
#                        experiment_target
#                    )
                if is_unique:
                    feat.regionID = reg_uid
                    feat.sampleID = sam_uid
                    feat.experimentID = exp.uid
                    feat.sourceID = sou.uid
#                    if feat_type == "histone":
#                        feat.histone_type = experiment_target
#                    if feat_type == "tf":
#                        feat.tf = experiment_target
                    session.add(feat)
                    session.commit()
#        # Do not cluster
#        else:
#            # For each accession, biosample...
#            for accession, biosample in metadata[(experiment_type, experiment_target)]:
#                # Load BED file
#                bed_obj = pybedtools.BedTool(
#                    os.path.join(exp_dummy_dir, "%s.bed" % accession))
#                # Get sample
#                sample = Sample()
#                if sample.is_unique(
#                    session,
#                    samples[biosample]["cell_or_tissue"],
#                    samples[biosample]["treatment"],
#                    samples[biosample]["cell_line"],
#                    samples[biosample]["cancer"]
#                ):
#                    sample.name = samples[biosample]["cell_or_tissue"]
#                    sample.treatment = samples[biosample]["treatment"]
#                    sample.cell_line = samples[biosample]["cell_line"]
#                    sample.cancer = samples[biosample]["cancer"]
#                    session.add(sample)
#                    session.commit()
#                # For each feature...
#                for feature in bed_obj:
#                    # Ignore non-standard chroms, scaffolds, etc.
#                    m = re.search("^chr(\S+)$", feature[0])
#                    if not m.group(1) in GUDglobals.chroms: continue
#                    # Get coordinates
#                    chrom = feature[0]
#                    start = int(feature[1])
#                    end = int(feature[2])
#                    # Get region
#                    region = Region()
#                    reg = region.select_by_exact_location(session, chrom, start, end)
#                    if not reg:
#                        # Insert region
#                        region.bin = assign_bin(start, end)
#                        region.chrom = chrom
#                        region.start = start
#                        region.end = end
#                        session.add(region)
#                        session.commit()
#                        reg = region.select_by_exact_location(session, chrom, start, end)
#                     # Insert feature
#                    if feat_type == "accessibility":
#                        feat = DNAAccessibility()
#                        is_unique = feat.is_unique(session,
#                            reg_uid, sou.uid, sam_uid, exp.uid)
#                    if feat_type == "histone":
#                        feat = HistoneModification()
#                        is_unique = feat.is_unique(session,
#                            reg_uid, sou.uid, sam_uid, exp.uid, experiment_target)
#                    if feat_type == "tf":
#                        feat = TFBinding()
#                        is_unique = feat.is_unique(session,
#                            reg_uid, sou.uid, sam_uid, exp.uid, experiment_target)
#                    if is_unique:
#                        feat.regionID = reg.uid
#                        feat.sourceID = sou.uid
#                        feat.sampleID = sam.uid
#                        feat.experimentID = exp.uid
#                        if feat_type == "histone":
#                            feat.histone_type = experiment_target
#                        if feat_type == "tf":
#                            feat.tf = experiment_target
#                        session.add(feat)
#                        session.commit()
#        # Remove dummy dir
#        if os.path.isdir(exp_dummy_dir): shutil.rmtree(exp_dummy_dir)

#-------------#
# Main        #
#-------------#

if __name__ == "__main__": main()