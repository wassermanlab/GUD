"""
Genomic Universal Database (GUD) module
"""

__author__ = "Oriol Fornes"
__credits__ = [
    "Oriol Fornes",
    "Tamar V. Av-Shalom",
    "Rachelle A. Farkas",
    "David J. Arenillas",
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
from sqlalchemy import create_engine
from sqlalchemy.orm import Session
import sys

__all__ = ["ORM", "scripts", "parsers"]

class Globals(object):
    """
    This class contains functions designed to
    work through the entire module.
    """

    #--------------#
    # Defaults     #
    #--------------#

    # Valid chromosomes
    chroms = list(map(str, range(1, 23))) + \
        ["X", "Y", "M"]

    # Defaults for the GUD database hosted at
    # the Wasserman lab
    db_name = "hg19"
    db_host = "ontarget.cmmt.ubc.ca"
    db_port = 5506
    db_user = "ontarget_r"

    # Defaults for gene2sample & sample2gene 
    max_samples = 0
    max_tss = 50
    min_percent = 0
    min_tpm = 100.0

    # Valid experiments
    experiments = [
        "ATAC-seq",
        "CAGE",
        "ChIP-seq",
        "DNase-seq",
        "FAIRE-seq",
        "GRO-seq"
    ]

    # Valid histone marks
    histone_types = [
        "H2A.Z",
        "H3K4me1",
        "H3K4me2",
        "H3K4me3",
        "H3K9ac",
        "H3K9me1",
        "H3K9me3",
        "H3K27ac",
        "H3K27me3",
        "H3K36me3",
        "H3K79me2",
        "H4K20me1"
    ]

    # Defaults for indentifying conserved
    # regions
    #
    # From Portales-Casamar et al. 2009:
    min_conservation = 0.7
    min_conservation_length = 20
    
    # Defaults for indentifying regulatory
    # regions
    #
    # From Blackwood et al. 1998: Transcriptional
    # control regions often contain multiple, au-
    # tonomous enhancer modules that vary from
    # about 50 bp to 1.5 kb in size
    #
    # From Sharan et al. 2007: Promoter length
    # can vary from 100 to 1000 base pairs
    min_enhancer_length = 100
    perc_score_thresh = 0.05

    # Valid repeat classes
    repeat_classes = ["LINE", "SINE", "LTR"]

    # Weights
    weights = {
        # General
        "conserved_regions": 1.0,
        "dna_accessibility": 1.0,
        "histone_modification": 1.0,
        "tf_binding": 1.0,
        # DNA accessibility
        "atac-seq": 1.0,
        "dnase-seq": 1.0,
        # Histone modifications
        "h2a.z": 0.0,
        "h3k4me1": 1.0,
        "h3k4me2": 1.0,
        "h3k4me3": 1.0,
        "h3k9ac": 1.0,
        "h3k9me1": 0.0,
        "h3k9me3": -1.0,
        "h3k27ac": 1.0,
        "h3k27me3": -1.0,
        "h3k36me3": 0.0,
        "h3k79me2": 0.0,
        "h4k20me1": 0.0
    }

    #--------------#
    # Input/Output #
    #--------------#

    def _get_file_handle(self,
        file_name, gz=False, mode="rt"):

        # Initialize
        raiseValueError = False
        
        # Open file handle
        if gz or file_name.endswith(".gz"):
            try:
                fh = gzip.open(file_name, mode)
            except:
                raiseValueError = True
        else:
            try:
                fh = open(file_name, mode)
            except:
                raiseValueError = True
        
        if raiseValueError:
            raise ValueError(
                "Could not open file handle: %s" % \
                file_name
            )

        return fh

    def parse_file(self, file_name, gz=False):
        """
        This function parses a file and yields lines
        one by one.

        @input:
        file_name {str}
        gz {bool} use the gzip module

        @yield: {str}
        """

        # Get file handle
        fh = self._get_file_handle(file_name, gz)

        # For each line...
        for line in fh:
            yield line.strip("\n")

        fh.close()

    def parse_csv_file(self, file_name, gz=False,
        delimiter=","):
        """
        This function parses a CSV file and yields
        lines one by one as a list.

        @input:
        file_name {str}
        gz {bool} use the gzip module
        delimiter {str} e.g. "\t"; default = ","

        @yield: {list}
        """

        # Get file handle
        fh = self._get_file_handle(file_name, gz)

        # For each line...
        for line in csv.reader(
            fh, delimiter=delimiter
        ):
            yield line

        fh.close()

    def parse_tsv_file(self, file_name, gz=False):
        """
        This function parses a tab-separated values
        file and yields lines one by one as a list.

        @input:
        file_name {str}
        gz {bool} use the gzip module

        @yield: {list}
        """

        # For each line...
        for line in self.parse_csv_file(
            file_name, gz, delimiter="\t"
        ):
            yield line

    def parse_fasta_file(self, file_name, gz=False):
        """
        This function parses a FASTA file and yields
        sequences one by one as a list of length 2
        (i.e. [{header}, {sequence}]).

        @input:
        file_name {str}
        gz {bool} use the gzip module to open the file

        @yield: {list}
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

    def write(self, file_name=None, content=None):
        """
        This function writes content to a file or,
        if no file is provided, to STDOUT. Content
        will be appended at the end of the file.
        """

        if file_name:
            # Get file handle
            fh = self._get_file_handle(
                file_name, mode="a"
            )
            # Write
            fh.write("%s\n" % content)
            fh.close()
        else:
            sys.stdout.write("%s\n" % content)

    #--------------#
    # SQLalchemy   #
    #--------------#

    def establish_GUD_session(self,
        db_user="ontarget_r",
        db_pass=None,
        db_host="ontarget.cmmt.ubc.ca",
        db_port=5506,
        db_name="hg19"
    ):

        if not db_pass: db_pass = ""
    
        gud_db = "mysql://{}:{}@{}:{}/{}".format(
            db_user,
            db_pass,
            db_host,
            db_port,
            db_name
        )

        # Establish a MySQL session
        try:
            engine = create_engine(
                gud_db,
                echo=False,
                pool_pre_ping=True
            )
            session = Session(engine)
        except:
            raise ValueError(
                "Could not connect to GUD db: %s" \
                % gud_db
            )

        return session

GUDglobals = Globals()