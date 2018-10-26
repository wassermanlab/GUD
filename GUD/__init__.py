"""
@file: __init__.py
@author: Oriol Fornes
@mail: oriol.fornes@gmail.com
@date: 2018
@ [wasserman lab](http://www.cisreg.ca)
@class: Globals
"""

import os, sys, re
import gzip

class Globals(object):
    """
    This class contains functions that have been designed to work
    through the whole {GUD} library.
    """

    #-------------#
    # Parsers     #
    #-------------#

    def parse_file(file_name, gz=False):
        """
        This function parses any file and yields lines one by one.

        @input:
        file_name {string}
        gz {boolean} if true, use the gzip module

        @return:
        line {string}
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

    def parse_csv_file(file_name, gz=False):
        """
        This function parses any CSV file and yields lines as a list.

        @input:
        file_name {string}
        gz {boolean} if true, use the gzip module

        @return: {list}
        line {list}
        """

        # For each line... #
        for line in parse_file(file_name, gz):
            line = line.split(",")
            yield line

    def parse_tsv_file(file_name, gz=False):
        """
        This function parses any TSV file and yields lines as a list.

        @input:
        file_name {string}
        gz {boolean} if true, use the gzip module

        @return: {list}
        line {list}
        """

        # For each line... #
        for line in parse_file(file_name, gz):
            line = line.split("\t")
            yield line

    def parse_fasta_file(file_name, gz=False, clean=True):
        """
        This function parses any FASTA file and yields sequences one
        by one as a list in the form [header, sequence].

        @input:
        file_name {string}
        gz {boolean} if true, use the gzip module
        clean {boolean} if true, convert non-amino acid res. to Xs

        @return:
        line {list} header, sequence
        """

        # Initialize #
        header = ""
        sequence = ""
        # For each line... #
        for line in parse_file(file_name, gz):
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

    #-------------#
    # Write       #
    #-------------#

    def write(file_name=None, content=None):
        """
        This function writes any {content} to a file or to stdout if no
        file is provided. If the file already exists, it pushed the {content}
        at the bottom of the file.

        @input:
        file_name {string}
        content {string}
        """
        if file_name is not None:
            try:
                f = open(file_name, "a")
                f.write("%s\n" % content)
            except:
                raise ValueError("Could create file %s" % file_name)
        else:
            sys.stdout.write("%s\n" % content)

GUDglobals = Globals()