"""
Genomic Universal Database (GUD) module
"""

__author__ = "Oriol Fornes"
__credits__ = ["Oriol Fornes", "Phillip A. Richmond", "David J. Arenillas",
    "Rachelle A. Farkas", "Michelle Kang", "Wyeth W. Wasserman"]
__email__ = "oriol@cmmt.ubc.ca"
__organization__ = "[Wasserman Lab](http://www.cisreg.ca)"
__version__ = "0.0.1"

import os, sys, re
import gzip

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

#    def parse_csv_file(self, file_name, gz=False):
#        """
#        This function parses a CSV file and yields lines one by one
#        as a list.
#
#        @input:
#        file_name {str}
#        gz {bool} use the gzip module
#
#        @return: {list}
#        """
#
#        # For each line... #
#        for line in self.parse_file(file_name, gz):
#            line = line.split(",")
#            yield line

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

#    def parse_fasta_file(self, file_name, gz=False, clean=True):
#        """
#        This function parses a FASTA file and yields sequences one
#        by one as a list in the form [header, sequence].
#
#        @input:
#        file_name {str}
#        gz {bool} use the gzip module
#
#        @return: [header, sequence]
#        """
#
#        # Initialize #
#        header = ""
#        sequence = ""
#        # For each line... #
#        for line in self.parse_file(file_name, gz):
#            if len(line) == 0: continue
#            if line.startswith("#"): continue
#            if line.startswith(">"):
#                if header != "" and sequence != "":
#                    yield header, sequence
#                header = ""
#                sequence = ""
#                m = re.search("^>(.+)", line)
#                if m: header = m.group(1)
#            elif header != "":
#                sub_sequence = line.upper()
#                if clean: sub_sequence = re.sub("[^ACDEFGHIKLMNPQRSTUVWY]", "X", sub_sequence)
#                sequence += sub_sequence
#        if header != "" and sequence != "":
#            yield header, sequence
#
#    #-------------#
#    # Write       #
#    #-------------#
#
#    def write(self, file_name=None, content=None):
#        """
#        This function writes {content} to a file or to stdout if no
#        file is provided. Note that {content} will be appended at
#        the end of the file.
#
#        @input:
#        file_name {str}
#        content {str}
#        """
#        if file_name is not None:
#            try:
#                f = open(file_name, "a")
#                f.write("%s\n" % content)
#            except:
#                raise ValueError("Could create file %s" % file_name)
#        else:
#            sys.stdout.write("%s\n" % content)

GUDglobals = Globals()
