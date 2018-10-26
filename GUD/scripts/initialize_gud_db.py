#!/usr/bin/env python2.7

import os, sys, re
import argparse
import ConfigParser
from datetime import date
from ftplib import FTP
import getpass
import gzip
from io import BytesIO
from sqlalchemy import create_engine
from sqlalchemy.orm import scoped_session, sessionmaker
from sqlalchemy_utils import create_database, database_exists

# Import from GUD module
from GUD.bin_range import BinRange
from GUD.ORM.chrom_size import ChromSize
from GUD.ORM.conservation import Conservation
from GUD.ORM.dna_accessibility import DnaAccessibility
from GUD.ORM.enhancer import Enhancer
from GUD.ORM.gene import Gene
from GUD.ORM.histone_modification import HistoneModification
from GUD.ORM.repeat_mask import RepeatMask
from GUD.ORM.tad import Tad
from GUD.ORM.tf_binding import TfBinding
#from GUD.ORM.tss import TSS

#-------------#
# Functions   #
#-------------#

def parse_args():
    """
    This function parses arguments provided via the command
    line using argparse.
    """

    parser = argparse.ArgumentParser(description="describe what the script does...")

    parser.add_argument("genome", help="Genome assembly")
    
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
    if args.db is None:
        args.db = args.genome

    return args

def initialize_gud_db(user, host, port, db, genome):
    """
    """

    # Initialize
    chroms = list(map(str, range(1, 23))) + ["X", "Y", "M"]
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
    
    # Create chrom sizes table
    if not engine.has_table("chrom_size"):
        # Initialize
        rows = []
        table = ChromSize()
        table.metadata.bind = engine
        # Create table
        table.metadata.create_all(engine)
        # Get UCSC FTP file
        directory, file_name = get_ftp_dir_and_file(genome, "chrom_size")
        # Download data
        for line in fetch_lines_from_ftp_file(
            genome, directory, file_name):
            # Split line
            line = line.split("\t")
            # Ignore non-standard chroms, scaffolds, etc.
            m = re.search("^chr(\S+)$", line[0])
            if not m.group(1) in chroms: continue
            # Add row
            rows.append(
                {
                    "chrom": line[0],
                    "size": line[1]
                }
            )
        # Insert rows to table
        engine.execute(table.__table__.insert(), rows)
    
    # Create conservation table
    if not engine.has_table("conservation"):
        # Initialize
        rows = []
        table = Conservation()
        table.metadata.bind = engine
        # Create table
        table.metadata.create_all(engine)
        # Get UCSC FTP file
        directory, file_name = get_ftp_dir_and_file(genome, "conservation")
        # Get source name
        source_name = re.search("^.+/(.+).txt.gz$", file_name)
        # Download data
        for line in fetch_lines_from_ftp_file(
            genome, directory, file_name):
            # Split line
            line = line.split("\t")
            # Ignore non-standard chroms, scaffolds, etc.
            m = re.search("^chr(\S+)$", line[1])
            if not m.group(1) in chroms: continue
            # Get bin
            start = int(line[2])
            end = int(line[3])
            bin = BinRange().binFromRange(start, end)
            # Add row
            rows.append(
                {
                    "bin": bin,
                    "chrom": line[1],
                    "chromStart": start,
                    "chromEnd": end,
#                    "extFile": line[4],
#                    "offset": line[5],
                    "score": line[6],
                    "source_name": source_name.group(1),
                    "date": today
                }
            )
            # Insert rows in bulks of 100,000
            if len(rows) == 100000:
                engine.execute(table.__table__.insert(), rows)
                # Clear rows
                rows = []
        # Insert remaining rows
        engine.execute(table.__table__.insert(), rows)

    # Create DNA accessibility table
    if not engine.has_table("dna_accessibility"):
        # Initialize
        rows = []
        table = DnaAccessibility()
        table.metadata.bind = engine
        # Create table
        table.metadata.create_all(engine)

    # Create enhancers table
    if not engine.has_table("enhancer"):
        # Initialize
        rows = []
        table = Enhancer()
        table.metadata.bind = engine
        # Create table
        table.metadata.create_all(engine)

    # Create gene table
    if not engine.has_table("gene"):
        # Initialize
        rows = []
        table = Gene()
        table.metadata.bind = engine
        # Create table
        table.metadata.create_all(engine)
        # Get UCSC FTP file
        directory, file_name = get_ftp_dir_and_file(genome, "gene")
        # Download data
        for line in fetch_lines_from_ftp_file(
            genome, directory, file_name):
            # Split line
            line = line.split("\t")
            # Ignore non-standard chroms, scaffolds, etc.
            m = re.search("^chr(\S+)$", line[2])
            if not m.group(1) in chroms: continue
            # Get bin
            start = int(line[2])
            end = int(line[3])
            bin = BinRange().binFromRange(start, end)
            # Add row
            rows.append(
                {
                    "bin": bin,
                    "name": line[1],
                    "chrom": line[2],
                    "strand": line[3],
                    "txStart": start,
                    "txEnd": end,
                    "cdsStart": line[6],
                    "cdsEnd": line[7],
#                    "exonCounts": line[8],
                    "exonStarts": line[9],
                    "exonEnds": line[10],
#                    "score": line[11],
                    "name2": line[12],
#                    "cdsStartStat": line[13],
#                    "cdsEndStat": line[14],
#                    "exonFrames": line[15],
                    "source_name": "refGene",
                    "date": today
                }
            )
        # Insert rows
        engine.execute(table.__table__.insert(), rows)

    # Create histone marks table
    if not engine.has_table("histone_modification"):
        # Initialize
        rows = []
        table = HistoneModification()
        table.metadata.bind = engine
        # Create table
        table.metadata.create_all(engine)

    # Create repeat mask table
    if not engine.has_table("rmsk"):
        # Initialize
        rows = []
        table = RepeatMask()
        table.metadata.bind = engine
        # Create table
        table.metadata.create_all(engine)
        # Get UCSC FTP file
        directory, file_name = get_ftp_dir_and_file(genome, "rmsk")
        # Download data
        for line in fetch_lines_from_ftp_file(
            genome, directory, file_name):
            # Split line
            line = line.split("\t")
            # Ignore non-standard chroms, scaffolds, etc.
            m = re.search("^chr(\S+)$", line[5])
            if not m.group(1) in chroms: continue
            # Get bin
            start = int(line[2])
            end = int(line[3])
            bin = BinRange().binFromRange(start, end)
            # Add row
            rows.append(
                {
                    "bin": bin,
                    "swScore": line[1],
#                    "milliDiv": line[2],
#                    "milliDel": line[3],
#                    "milliIns": line[4],
                    "genoName": line[5],
                    "genoStart": start,
                    "genoEnd": end,
#                    "genoLeft": line[8],
                    "strand": line[9],
                    "repName": line[10],
                    "repClass": line[11],
                    "repFamily": line[12],
#                    "repStart": line[13],
#                    "repEnd": line[14],
#                    "repLeft": line[15],
#                    "id": line[16],
                    "source_name": "rmsk",
                    "date": today
                }
            )
            # Insert rows in bulks of 100,000
            if len(rows) == 100000:
                engine.execute(table.__table__.insert(), rows)
                # Clear rows
                rows = []
        # Insert remaining rows
        engine.execute(table.__table__.insert(), rows)            

    # Create TAD table
    if not engine.has_table("tad"):
        # Initialize
        rows = []
        table = Tad()
        table.metadata.bind = engine
        # Create table
        table.metadata.create_all(engine)

    # Create TF-binding table
    if not engine.has_table("tf_binding"):
        # Initialize
        rows = []
        table = TfBinding()
        table.metadata.bind = engine
        # Create table
        table.metadata.create_all(engine)

#    # Create TSS table
#    if not engine.has_table("tss"):
#        # TO BE DONE!!!
#        # i.e. the FANTOM5 table!

def get_ftp_dir_and_file(genome, data_type):

    # Initialize
    ftp = FTP("hgdownload.soe.ucsc.edu")
    ftp.login()

    # Change into "genome" directory
    try:
        ftp.cwd(os.path.join("goldenPath", genome))
    except:
        raise ValueError("Cannot connect to FTP goldenPath site: %s" % genome)

    # Fetch bigZips and database files
    if data_type == "chrom_size":
        return "bigZips", "%s.chrom.sizes" % genome
    elif data_type == "gene":
        return "database", "refGene.txt.gz"
    elif data_type == "rmsk":
        return "database", "rmsk.txt.gz"
    elif data_type == "conservation":
        regexp = re.compile("multiz\d+way.txt.gz$")
        for file_name in sorted(filter(regexp.search, ftp.nlst("database"))):
            return "database", file_name

def fetch_lines_from_ftp_file(genome, directory, file_name):
    
    # Initialize
    global BIO
    ftp = FTP("hgdownload.soe.ucsc.edu")
    ftp.login()
    BIO = BytesIO()

    # Change into "genome" directory
    try:
        ftp.cwd(os.path.join("goldenPath", genome, directory))
    except:
        raise ValueError("Cannot connect to FTP goldenPath site: %s/%s" % (genome, directory))

    # If valid file...
    if file_name in ftp.nlst():
        # Retrieve FTP file
        ftp.retrbinary("RETR %s" % file_name, callback=handle_bytes)
        BIO.seek(0) # Go back to the start
        # If compressed file...
        if file_name.endswith(".gz"):
            f = gzip.GzipFile(fileobj=BIO, mode="rb")
        # ... Else...
        else:
            f = BIO
        # For each line...
        for line in f:
            yield line.decode("UTF-8").strip("\n")

def handle_bytes(bytes):
    BIO.write(bytes)

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    # Parse arguments
    args = parse_args()

    # Initialize GUD database: create tables
    # and download data to populate them 
    initialize_gud_db(args.user, args.host,
        args.port, args.db, args.genome)
