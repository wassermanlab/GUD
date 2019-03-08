"""
Genomic Universal Database (GUD) module
"""

__author__ = "Oriol Fornes"
__credits__ = [
    "Oriol Fornes",
    "Tamar V. Av-Shalom",
    "David J. Arenillas",
    "Rachelle A. Farkas",
    "Michelle Kang",
    "Phillip A. Richmond",
    "Wyeth W. Wasserman"
]
__email__ = "oriol@cmmt.ubc.ca"
__organization__ = "[Wasserman Lab](http://www.cisreg.ca)"
__version__ = "0.0.1"

import os, sys, re
from Bio import SeqIO
import csv
import gzip

__all__ = ["ORM", "scripts"]

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
            # Open file handle
            if gz or file_name.endswith(".gz"):
                try: f = gzip.open(file_name, "r")
                except: raise ValueError("Could not open file %s" % file_name)
            else:
                try: f = open(file_name, "r")
                except: raise ValueError("Could not open file %s" % file_name)
            # For each line...
            for line in f:
                yield line.replace("\r", "").strip("\n")
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

        # For each line...
        for line in self.parse_file(file_name, gz):
            yield line.split(",")

    def parse_tsv_file(self, file_name, gz=False):
        """
        This function parses a TSV file and yields lines one by one
        as a list.

        @input:
        file_name {str}
        gz {bool} use the gzip module

        @return: {list}
        """

        # For each line...
        for line in self.parse_file(file_name, gz):
            yield line.split("\t")

    def parse_fasta_file(self, file_name, gz=False, clean=True):
        """
        This function parses a FASTA file and yields sequences one
        by one as a list in the form [header, sequence].

        @input:
        file_name {str}
        gz {bool} use the gzip module
        clean {bool} replace non-standard amino acids by Xs

        @return: [header, sequence]
        """

        if os.path.exists(file_name):
            # Open file handle
            if gz or file_name.endswith(".gz"):
                try: f = gzip.open(file_name, "r")
                except: raise ValueError("Could not open file %s" % file_name)
            else:
                try: f = open(file_name, "r")
                except: raise ValueError("Could not open file %s" % file_name)
            # For each SeqRecord...
            for seq_record in SeqIO.parse(f, "fasta"):
                # Initialize
                header = seq_record.id
                sequence = str(seq_record.seq).upper()
                yield header, sequence
            f.close()
        else:
            raise ValueError("File %s does not exist!" % file_name)

    #-------------#
    # Functions   #
    #-------------#

    def write(self, file_name=None, content=None):
        """
        This function writes content to a file or, if no file is
        provided, to STDOUT. Note that content will be appended
        at the end of the file.
        """

        if file_name:
            try:
                f = open(file_name, "a")
            except:
                raise ValueError("Could not write to file: %s" % file_name)
            # Write
            f.write("%s\n" % content)
            f.close()
        else:
            sys.stdout.write("%s\n" % content)

GUDglobals = Globals()