#!/usr/bin/env python
# mysql -u ontarget_w --database tamar_test

import re
import os
import sys
import re
import argparse
from binning import assign_bin
import getpass
from sqlalchemy import create_engine
from sqlalchemy.orm import mapper, scoped_session, sessionmaker
from sqlalchemy_utils import database_exists
import warnings
from datetime import date

# Import from GUD module
from GUD2 import GUDglobals
from GUD2.ORM.clinvar import ClinVar
from GUD2.ORM.region import Region
from GUD2.ORM.source import Source

#-------------#
# Functions   #
#-------------#
def parse_args():
    """
    This function parses arguments provided via the command
    line using argparse.
    """

    parser = argparse.ArgumentParser(
        description="this script inserts clinvar information")

    parser.add_argument("vcf_file", help="annotated clinvar")

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


def insert_clinvar_to_gud_db(user, host, port, db, vcf_file):
    # Initialize
    metadata = {}
    db_name = "mysql://{}:@{}:{}/{}".format(
        user, host, port, db)
    if not database_exists(db_name):
        raise ValueError("GUD db does not exist: %s" % db_name)
    session = scoped_session(sessionmaker())
    engine = create_engine(db_name, echo=False)
    session.remove()
    session.configure(bind=engine, autoflush=True,
                      expire_on_commit=False)

    # Initialize table
    table = ClinVar()
    table.metadata.bind = engine
    try:
        table.metadata.create_all(engine)
    except:
        raise ValueError("Cannot create \"ClinVar\" table!")
    
    # add source  
    day = str(date.today())
    with open(vcf_file) as f:
        for line in f : 
          if line.startswith("##fileDate="):
              day = line.split("=")[1].rstrip()
    source_name = "ClinVar_" + str(day)
    source = Source()
    sou = source.select_by_name(session, source_name)
    if not sou: 
        source.name = source_name
        session.merge(source)
        session.commit()
        sou = source.select_by_name(session, source_name)

    #parse table 
    annotation_list = ["ANN_Allele", "ANN_Annotation", 
    "ANN_Annotation_Impact", "ANN_Gene_Name", "ANN_Gene_ID", 
    "ANN_Feature_Type", "ANN_Feature_ID"]
    info_list = ["CADD","CLNDISDB","CLNDN","CLNSIG",
    "gnomad_exome_af_global","gnomad_exome_hom_global",
    "gnomad_genome_af_global","gnomad_genome_hom_global"]
    columns = annotation_list + info_list 
    with open(vcf_file) as f:
      for line in f:
        if not line.startswith("#"):
          line_list = [None]*len(columns)
          fields = line.split("\t")
          info = fields[-1].rstrip()
          info = info.split(";")
          info[-1] = info[-1].rstrip()
          
          # make region
          chrom = fields[0]
          start = int(fields[1]) - 1
          end = start + len(fields[3]) ##finish this 
          region = Region()
          reg = region.select_by_exact_location(session, chrom, start, end)
          if not reg: 
              region.bin = assign_bin(start, end)
              region.chrom = chrom
              region.start = start
              region.end = end 
              session.merge(region)
              session.commit()
              reg = region.select_by_pos(session, chrom, start, end)
          ## add info fields
          for i in info:
            i_split = i.split("=")
            if len(i_split) == 2:
              if i_split[0] in info_list:
                index = columns.index(i_split[0])
                line_list[index] = i_split[1].rstrip()
              if i_split[0] == "ANN": ## further split the annotation column
                ann = i_split[1].split(",")[0].split("|")
                for idx, val in enumerate(annotation_list):
                  index = columns.index(val)
                  line_list[index] = ann[idx]
          ## make clinvar 
          clinvarID = fields[2]
          clinvar = ClinVar()
          cln = clinvar.is_unique(session, clinvarID)
          if not cln: 
            clinvar.regionID = reg.uid
            clinvar.sourceID = sou.uid
            clinvar.ref = fields[3]
            clinvar.alt = fields[4]
            clinvar.clinvarID = fields[2]
            clinvar.ANN_Annotation = line_list[1]
            clinvar.ANN_Annotation_Impact = line_list[2]
            clinvar.ANN_Gene_Name = line_list[3]
            clinvar.ANN_Gene_ID = line_list[4]
            clinvar.ANN_Feature_Type = line_list[5]
            clinvar.ANN_Feature_ID = line_list[6]
            clinvar.CADD = None if line_list[7] is None else float(line_list[7])
            clinvar.CLNDISDB = line_list[8]
            clinvar.CLNDN = line_list[9]
            clinvar.CLNSIG = line_list[10]
            clinvar.gnomad_exome_af_global = None if line_list[11] is None else float(line_list[11])
            clinvar.gnomad_exome_hom_global = None if line_list[12] is None else float(line_list[12])
            clinvar.gnomad_genome_af_global = None if line_list[13] is None else float(line_list[13])
            clinvar.gnomad_genome_hom_global = None if line_list[14] is None else float(line_list[14])
          
#-------------#
# Main        #
#-------------#


if __name__ == "__main__":

    # Parse arguments
    args = parse_args()
    # Insert ENCODE data to GUD database
    insert_clinvar_to_gud_db(args.user, args.host, args.port, args.db,
                         args.vcf_file)