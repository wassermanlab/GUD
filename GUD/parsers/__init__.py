"""
Contains Python parsers with which to populate GUD
"""

from Bio import SeqIO
from io import BytesIO
from ftplib import FTP
import gzip
import os
import pandas
import re
from sqlalchemy_utils import create_database, database_exists
import sys
from zipfile import ZipFile

from GUD.ORM.chrom import Chrom
from GUD.ORM.conservation import Conservation
from GUD.ORM.dna_accessibility import DNAAccessibility
from GUD.ORM.experiment import Experiment
from GUD.ORM.gene import Gene
from GUD.ORM.histone_modification import HistoneModification
from GUD.ORM.mask import Mask
from GUD.ORM.region import Region
from GUD.ORM.sample import Sample
from GUD.ORM.source import Source
from GUD.ORM.tf_binding import TFBinding

class ParseUtililities:
    """
    Contains functions designed to work with the parsers.
    """

    #--------------#
    # Defaults     #
    #--------------#

    def __init__(self):
        """
        @param genome = genome assembly
        @param dbname = database name (i.e. from GUDUtils._get_db_name())
        @param engine = SQLAlchemy {Engine}
        """

        self._genome = None
        self._dbname = None
        self._engine = None
        self._session = None

    @property
    def genome(self):
        """
        @rtype = {String}
        """
        return(self._genome)

    @genome.setter
    def genome(self, value):
        self._genome = str(value)

    @property
    def dbname(self):
        """
        @rtype = {String}
        """
        return(self._dbname)

    @dbname.setter
    def dbname(self, value):
        self._dbname = str(value)

    # Engine-specific
    @property
    def engine(self):
        """
        From https://docs.sqlalchemy.org/en/13/core/connections.html
        @rtype = {Engine}
        """
        return(self._engine)

    @engine.setter
    def engine(self, value):
        self._engine = value

    def dispose_engine(self):
        self._engine.dispose()

    #--------------#
    # Input/Output #
    #--------------#

    def _get_file_handle(self, file_name, mode="r"):

        # Initialize
        raiseValueError = False
        
        # Open file handle
        if file_name.endswith(".gz"):
            try:
                handle = gzip.open(file_name, mode)
            except:
                raiseValueError = True

        elif file_name.endswith(".zip"):
            try:
                zf = ZipFile(file_name, mode)
                for f in zf.infolist():

                    handle = zf.open(f, mode)
                    break
            except:
                raiseValueError = True

        else:
            try:
                handle = open(file_name, mode)
            except:
                raiseValueError = True
        
        if raiseValueError:
            raise ValueError("Could not open file handle: %s" % file_name)

        return(handle)

    def parse_file(self, file_name):
        """
        Parses a file and yields lines one by one.

        @input:
        file_name {str}

        @yield: {str}
        """

        # Get file handle
        handle = self._get_file_handle(file_name)

        # For each line...
        for line in handle:

            yield(line.strip("\n"))

        handle.close()

    def parse_csv_file(self, file_name, delimiter=","):
        """
        Parses a CSV file and yields lines one by one as a list.

        @input:
        file_name {str}
        delimiter {str} e.g. "\t"; default = ","

        @yield: {list}
        """

        # Get file handle
        handle = self._get_file_handle(file_name)

        # Read in chunks
        for chunk in pandas.read_csv(handle, header=None, encoding="utf8", sep=delimiter, chunksize=1024):
            for index, row in chunk.iterrows():
                yield(row.tolist())

        handle.close()

    def parse_tsv_file(self, file_name):
        """
        Parses a TSV file and yields lines one by one as a list.

        @input:
        file_name {str}

        @yield: {list}
        """

        # For each line...
        for line in self.parse_csv_file(file_name, delimiter="\t"):

            yield(line)

    def parse_fasta_file(self, file_name):
        """
        Parses a FASTA file and yields sequences one by one  as a list of
        length 2 (i.e. [{header}, {sequence}]).

        @input:
        file_name {str}

        @yield: {list}
        """

        # Get file handle
        handle = self._get_file_handle(file_name)

        # For each SeqRecord...
        for seq_record in SeqIO.parse(handle, "fasta"):
            # Initialize
            header = seq_record.id
            sequence = str(seq_record.seq).upper()

            yield(header, sequence)

        handle.close()

    # def parse_file_in_chunks(self, file_name=None, chunk_size=1024):
    #     """
    #     Parses a file and yields lines in chunks.

    #     @input:
    #     file_name {str}
    #     chunk_size {int}

    #     @yield: {str}
    #     """

    #     # Get file handle
    #     handle = self._get_file_handle(file_name)

    #     while True:

    #         data = handle.read(chunk_size)

    #         if not data:
    #             break

    #         yield(data)

    #     handle.close()

    def write(self, file_name=None, content=None):
        """
        Writes content to a file or, if no file is provided, to STDOUT.
        Content will be appended at the end of the file.
        """

        if file_name:

            # Get file handle
            handle = self._get_file_handle(file_name, mode="a")

            # Write
            handle.write("%s\n" % content)
            handle.close()

        else:

            sys.stdout.write("%s\n" % content)

    #--------------#
    # Initializer  #
    #--------------#

    def initialize_gud_db(self):

        if not database_exists(self.dbname):

            # Initialize
            chroms = []
            valid_chroms = set(list(map(str, range(1, 23))) + ["X", "Y", "M"])

            # Create database
            create_database(self.dbname)

            # For each table...
            for table in [Chrom(), Experiment(), Region(), Sample(), Source()]:
                    self.create_table(table)
    
            # Download data
            for line in self._download_data_from_UCSC():

                # Split line
                line = line.split("\t")

                # Ignore non-standard chroms, scaffolds, etc.
                m = re.search("^chr(.+)$", line[0])
                if not m.group(1) in valid_chroms:
                    continue

                # Append chromosome
                chroms.append({"chrom": line[0], "size": line[1]})

            # Insert features
            self.engine.execute(Chrom().__table__.insert(), chroms)

    def create_table(self, table):

        if not self.engine.has_table(table.__tablename__):
            table.__table__.create(bind=self.engine)

    def _download_data_from_UCSC(self):

        # Initialize
        global BIO
        BIO = BytesIO()
        ftp_file = "%s.chrom.sizes" % self.genome

        # Login to UCSC's FTP server
        ftp = FTP("hgdownload.soe.ucsc.edu")
        ftp.login()

        # Change into bigZips folder for genome
        try:
            ftp.cwd(os.path.join("goldenPath", self.genome, "bigZips"))
        except:
            error = ["error", "cannot connect to UCSC's FTP goldenPath folder", self.genome, "bigZips"]
            print(": ".join(error))
            exit(0)

        # Retrieve FTP file
        ftp.retrbinary("RETR %s" % ftp_file, callback=self._handle_bytes)
        BIO.seek(0) # Go back to the start

        # For each line...
        for line in BIO:
            yield(line.decode("utf-8").strip("\n"))

    def _handle_bytes(self, bytes):
        BIO.write(bytes)

    #--------------#
    # Getters      #
    #--------------#

    def get_chroms(self, session):
        return(Chrom.chrom_sizes(session))

    def get_experiment(self, session, experiment_type):
        return(Experiment.select_unique(session, experiment_type))

    def get_region(self, session, chrom, start, end, strand=None):
        return(Region.select_unique(session, chrom, start, end, strand))

    def get_sample(self, session, name, X, Y, treatment, cell_line, cancer):
        """
        @classmethod
        def select_unique(cls, session, name, treatment, cell_line, cancer):
        """
        return(Sample.select_unique(session, name, X, Y, treatment, cell_line, cancer))

    def get_source(self, session, source_name):
        return(Source.select_unique(session, source_name))

    #--------------#
    # Upserts      #
    #--------------#

    def upsert_accessibility(self, session, accessibility):

        if DNAAccessibility.is_unique(session, accessibility.region_id, accessibility.sample_id, accessibility.experiment_id, accessibility.source_id):
            session.add(accessibility)
            session.flush()

    def upsert_conservation(self, session, conservation):

        if Conservation.is_unique(session, conservation.region_id, conservation.source_id):
            session.add(conservation)
            session.flush()

    def upsert_experiment(self, session, experiment):

        if Experiment.is_unique(session, experiment.name):
            session.add(experiment)
            session.flush()

    def upsert_gene(self, session, gene):

        if Gene.is_unique(session, gene.region_id, gene.name, gene.source_id):
            session.add(gene)
            session.flush()

    def upsert_histone(self, session, histone):

        if HistoneModification.is_unique(session, histone.region_id, histone.sample_id, histone.experiment_id, histone.source_id, histone.histone_type):
            session.add(histone)
            session.flush()

    def upsert_mask(self, session, mask):

        if Mask.is_unique(session, mask.region_id, mask.source_id):
            session.add(mask)
            session.flush()

    def upsert_region(self, session, region):

        if Region.is_unique(session, region.chrom, region.start, region.end, region.strand):
            session.add(region)
            session.flush()

    def upsert_sample(self, session, sample):
        """
        @classmethod
        def is_unique(cls, session, name, X, Y, treatment, cell_line, cancer):
        """
        if Sample.is_unique(session, sample.name, sample.X, sample.Y, sample.treatment, sample.cell_line, sample.cancer):
            session.add(sample)
            session.flush()

    def upsert_source(self, session, source):

        if Source.is_unique(session, source.name):
            session.add(source)
            session.flush()

    def upsert_tf(self, session, tf):

        if TFBinding.is_unique(session, tf.region_id, tf.sample_id, tf.experiment_id, tf.source_id, tf.tf):
            session.add(tf)
            session.flush()

    #--------------#
    # Multiprocess #
    #--------------#

    def insert_data_files_in_parallel(self, data_files, insert_function, threads=1):

        # from itertools import islice
        from multiprocessing import Pool

        # Parallelize
        pool = Pool(processes=threads)
        pool.map(insert_function, data_files)
        pool.close()

    #     # Get iterable
    #     iterable = self.parse_tsv_file(data_file)

    #     # Get chunks
    #     chunks = self._grouper(iterable)

    #     while True:

    #         # Groups chunks for multiprocessing
    #         grouped_chunks = [list(chunk) for chunk in islice(chunks, threads)]

    #         if grouped_chunks:
    #             pool.map(insert_function, grouped_chunks)

    #         else:
    #             break

    #         if test:
    #             break

    #     if test:
    #         for chunk in islice(chunks, threads):
    #             continue

    #     # Close pool
    #     pool.close()

    # def _grouper(self, iterable, n=5000, fillvalue=None):

    #     import sys
    #     # Python 3+
    #     if sys.version_info > (3, 0):
    #         from itertools import zip_longest
    #     # Python 2.7
    #     else:
    #         from itertools import izip_longest as zip_longest

    #     args = [iter(iterable)] * n

    #     return(zip_longest(*args, fillvalue=fillvalue))

ParseUtils = ParseUtililities()