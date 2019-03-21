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

from Bio import SeqIO
import csv
import gzip
import os
import sys

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

    def _get_file_handle(self, file_name, gz=False, mode="r"):

        if os.path.exists(file_name):
            # Open file handle
            if gz or file_name.endswith(".gz"):
                try:
                    fh = gzip.open(file_name, mode)
                except:
                    raise ValueError("Could not open file handle %s" % file_name)
            else:
                try:
                    fh = open(file_name, mode)
                except:
                    raise ValueError("Could not open file handle %s" % file_name)
        else:
            # File does not exist
            raise ValueError("Could not open file handle %s" % file_name)

        return fh

    def parse_file(self, file_name, gz=False):
        """
        This function parses a file and yields lines one by one.

        @input:
        file_name {str}
        gz {bool} use the gzip module to open the file

        @yield: {str}
        """

        # Get file handle
        fh = self._get_file_handle(file_name, gz)

        # For each line...
        for line in fh:
            yield line.replace("\r", "").strip("\n")

        fh.close()

    def parse_csv_file(self, file_name, gz=False, delimiter=","):
        """
        This function parses a CSV file and yields lines one by
        one as a list.

        @input:
        file_name {str}
        gz {bool} use the gzip module to open the file
        delimiter {str} e.g. "\t"; default = ","

        @yield: {list}
        """

        # Get file handle
        fh = self._get_file_handle(file_name, gz)

        # For each line...
        for line in csv.reader(fh, delimiter):
            yield line

        fh.close()

    def parse_tsv_file(self, file_name, gz=False):
        """
        This function parses a tab-separated values file and
        yields lines one by one as a list.

        @input:
        file_name {str}
        gz {bool} use the gzip module to open the file

        @yield: {list}
        """

        # For each line...
        for line in self.parse_csv_file(file_name, gz, delimiter="\t"):
            yield line

    def parse_fasta_file(self, file_name, gz=False, clean=True):
        """
        This function parses a FASTA file and yields sequences
        one by one as a list in the form [header, sequence].

        @input:
        file_name {str}
        gz {bool} use the gzip module to open the file
        clean {bool} replace non-standard amino acids by Xs

        @yield: [header, sequence]
        """

        # Get file handle
        fh = self._get_file_handle(file_name, gz)

        # For each SeqRecord...
            for seq_record in SeqIO.parse(fh, "fasta"):
                # Initialize
                header = seq_record.id
                sequence = str(seq_record.seq).upper()

                yield header, sequence

        fh.close()

    #-------------#
    # Functions   #
    #-------------#

    def write(self, file_name=None, content=None):
        """
        This function writes content to a file or, if no file
        file is provided, to STDOUT. Note that content will be
        appended at the end of the file.
        """

        if file_name:
            # Get file handle
            fh = self._get_file_handle(file_name, mode="a")
            # Write
            fh.write("%s\n" % content)
            fh.close()
        else:
            sys.stdout.write("%s\n" % content)

GUDglobals = Globals()