"""
Genomic Universal Database (GUD) module
"""

__author__ = "Oriol Fornes"
__credits__ = ["Oriol Fornes", "Phillip A. Richmond", "Tamar V. Av-Shalom",
    "David J. Arenillas", "Rachelle A. Farkas", "Michelle Kang", 
    "Wyeth W. Wasserman"]
__email__ = "oriol@cmmt.ubc.ca"
__organization__ = "[Wasserman Lab](http://www.cisreg.ca)"
__version__ = "0.0.1"

import os, sys, re
from ftplib import FTP
import gzip
from io import BytesIO

__all__ = ["ORM"]

class Globals(object):
    """
    This class contains functions designed to work through the
    entire module.
    """

    #-------------#
    # Definitions #
    #-------------#

    chroms = list(map(str, range(1, 23))) + ["X", "Y", "M"]

    #-------------#
    # Parsers     #
    #-------------#

    def parse_file(self, file_name, gz=False):
        """
        This function parses a file and yields lines one by one.

        @input:
        file_name {str}
        gz {bool} use the gzip module

        @return: {str}
        """

        if os.path.exists(file_name):
            # Initialize #
            f = None
            # Open file handle #
            if gz:
                try: f = gzip.open(file_name, "rt")
                except: raise ValueError("Could not open file %s" % file_name)
            else:
                try: f = open(file_name, "rt")
                except: raise ValueError("Could not open file %s" % file_name)
            # For each line... #
            for line in f:
                line = line.replace("\r", "")
                yield line.strip("\n")
            f.close()
        else:
            raise ValueError("File %s does not exist!" % file_name)

    def parse_csv_file(self, file_name, gz=False):
        """
        This function parses a CSV file and yields lines one by one
        as a list.

        @input:
        file_name {str}
        gz {bool} use the gzip module

        @return: {list}
        """

        # For each line... #
        for line in self.parse_file(file_name, gz):
            line = line.split(",")
            yield line

    def parse_tsv_file(self, file_name, gz=False):
        """
        This function parses a TSV file and yields lines one by one
        as a list.

        @input:
        file_name {str}
        gz {bool} use the gzip module

        @return: {list}
        """

        # For each line... #
        for line in self.parse_file(file_name, gz):
            line = line.split("\t")
            yield line

    def parse_fasta_file(self, file_name, gz=False, clean=True):
        """
        This function parses a FASTA file and yields sequences one
        by one as a list in the form [header, sequence].

        @input:
        file_name {str}
        gz {bool} use the gzip module

        @return: [header, sequence]
        """

        # Initialize #
        header = ""
        sequence = ""
        # For each line... #
        for line in self.parse_file(file_name, gz):
            if len(line) == 0: continue
            if line.startswith("#"): continue
            if line.startswith(">"):
                if header != "" and sequence != "":
                    yield header, sequence
                header = ""
                sequence = ""
                m = re.search("^>(.+)", line)
                if m: header = m.group(1)
            elif header != "":
                sub_sequence = line.upper()
                if clean: sub_sequence = re.sub("[^ACDEFGHIKLMNPQRSTUVWY]", "X", sub_sequence)
                sequence += sub_sequence
        if header != "" and sequence != "":
            yield header, sequence

    def get_ucsc_ftp_dir_and_file(self, genome, data_type):

        # Initialize
        ftp = FTP("hgdownload.soe.ucsc.edu")
        ftp.login()

        # Change into "genome" folder
        try:
            ftp.cwd(os.path.join("goldenPath", genome))
        except:
            raise ValueError("Cannot connect to FTP goldenPath folder: %s" % genome)

        # Fetch bigZips and database files
        if data_type == "chrom_size":
            return "bigZips", "%s.chrom.sizes" % genome
        elif data_type == "gene":
            return "database", "refGene.txt.gz"
        elif data_type == "rmsk":
            return "database", "rmsk.txt.gz"
        elif data_type == "conservation":
            regexp = re.compile("(multiz\d+way.txt.gz)")
            for file_name in sorted(filter(regexp.search, ftp.nlst("database"))):
                m = re.search(regexp, file_name)
                return "database", m.group(1)

    def fetch_lines_from_ucsc_ftp_file(self, genome, directory, file_name):

        # Initialize
        global BIO
        ftp = FTP("hgdownload.soe.ucsc.edu")
        ftp.login()
        BIO = BytesIO()

        # Change into "genome" "directory" folder
        try:
            ftp.cwd(os.path.join("goldenPath", genome, directory))
        except:
            raise ValueError("Cannot connect to FTP goldenPath folder: %s/%s" % (genome, directory))

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

    def handle_bytes(self, bytes):
        BIO.write(bytes)

GUDglobals = Globals()
