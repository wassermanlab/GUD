"""
Contains Python parsers with which to populate GUD
"""

from io import BytesIO
from ftplib import FTP
import gzip
import os
import re
from sqlalchemy import create_engine
from sqlalchemy.orm import scoped_session, sessionmaker
from sqlalchemy_utils import create_database

from GUD import GUDglobals
from GUD.ORM.chrom import Chrom
from GUD.ORM.conservation import Conservation
from GUD.ORM.experiment import Experiment
from GUD.ORM.gene import Gene
from GUD.ORM.region import Region
from GUD.ORM.sample import Sample
from GUD.ORM.source import Source
from GUD.ORM.tf_binding import TFBinding

#--------------#
# Initializer  #
#--------------#

def _initialize_gud_db(user, pwd, host, port, db, genome):

    # Create database
    db_name = _get_db_name(user, pwd, host, port, db)
    create_database(db_name)

    # Get engine/session
    engine, Session = _initialize_engine_session(db_name)

    # For each table...
    for table in [Chrom(), Experiment(), Region(), Sample(), Source()]:

        if not engine.has_table(table.__tablename__):
            # Create table
            table.__table__.create(bind=engine)

            # If chromosomes table...
            if table.__tablename__ == "chroms":

                # Insert chrom sizes
                _insert_chrom_sizes(engine, table, genome)

    # Close connection
    Session.remove()
    engine.dispose()

def _insert_chrom_sizes(engine, table, genome):

    # Intialize
    chroms = []

    # Download data
    for line in _download_data_from_UCSC(genome):

        # Split line
        line = line.split("\t")

        # Ignore non-standard chroms, scaffolds, etc.
        m = re.search("^chr(.+)$", line[0])
        if not m.group(1) in GUDglobals.chroms:
            continue

        # Append chromosome
        chroms.append({"chrom": line[0], "size": line[1]})

    # Insert features
    engine.execute(table.__table__.insert(), chroms)

def _download_data_from_UCSC(genome):

    # Initialize
    global BIO
    BIO = BytesIO()
    ftp_file = "%s.chrom.sizes" % genome

    # Login to UCSC's FTP server
    ftp = FTP("hgdownload.soe.ucsc.edu")
    ftp.login()

    # Change into bigZips folder for genome
    try:
        ftp.cwd(os.path.join("goldenPath", genome, "bigZips"))
    except:
        error = ["error", "cannot connect to UCSC's FTP goldenPath folder", genome, "bigZips"]
        print(": ".join(error))
        exit(0)

    # Retrieve FTP file
    ftp.retrbinary("RETR %s" % ftp_file, callback=_handle_bytes)
    BIO.seek(0) # Go back to the start

    # For each line...
    for line in BIO:
        yield(line.decode("utf-8").strip("\n"))

    # # Fetch bigZips and database files
    # if data_type == "chroms":
    #     directory = "bigZips"
    #     file_name = 
    # if data_type == "conservation":
    #     regexp = re.compile("(multiz[0-9]+way.txt.gz)")
    #     for file_name in sorted(filter(regexp.search, ftp.nlst("database"))):
    #         m = re.search(regexp, file_name)
    #         directory = "database"
    #         file_name = m.group(1)
    #         break
    # if data_type == "gene":
    #     directory = "database"
    #     file_name = "refGene.txt.gz"

    # # Change into "directory" folder
    # try:
    #     ftp.cwd(directory)
    # except:
    #     error = ["error", "cannot connect to UCSC's FTP goldenPath folder", genome, directory]
    #     print(": ".join(error))
    #     exit(0)

    # # If file exists...
    # if file_name in ftp.nlst():

    #     # Retrieve FTP file
    #     ftp.retrbinary("RETR %s" % file_name, callback=_handle_bytes)
    #     BIO.seek(0) # Go back to the start

    #     # If compressed file...
    #     if file_name.endswith(".gz"):
    #         f = gzip.GzipFile(fileobj=BIO, mode="rb")
    #     else:
    #         f = BIO

    #     # For each line...
    #     for line in f:
    #         yield(line.decode("utf-8").strip("\n"))

def _handle_bytes(bytes):
    BIO.write(bytes)


#--------------#
# SqlAlchemy   #
#--------------#

def _get_db_name(user, pwd, host, port, db):
    return("mysql+pymysql://{}:{}@{}:{}/{}".format(user, pwd, host, port, db))

def _initialize_engine_session(db_name):

    # Initialize
    engine = create_engine(db_name, pool_pre_ping=True, pool_size=20, max_overflow=0)
    session_factory = sessionmaker(bind=engine)
    Session = scoped_session(session_factory)

    return(engine, Session)

#--------------#
# Getters      #
#--------------#

def _get_chroms(session):
    return(Chrom.chrom_sizes(session))

def _get_experiment(session, experiment_type):
    return(Experiment.select_by_name(session, experiment_type))

def _get_region(session, chrom, start, end, strand=None):
    return(Region.select_unique(session, chrom, start, end, strand))

def _get_sample(session, name, treatment, cell_line, cancer):
    """
    @classmethod
    def select_unique(cls, session, name, treatment, cell_line, cancer):
    """
    return(Sample.select_unique(session, name, treatment, cell_line, cancer))

def _get_source(session, source_name):
    return(Source.select_by_name(session, source_name))

#--------------#
# Upserts      #
#--------------#

def _upsert_conservation(session, conservation):

    if Conservation.is_unique(session, conservation.regionID, conservation.sourceID):
        session.add(conservation)
        session.flush()

def _upsert_experiment(session, experiment):

    if Experiment.is_unique(session, experiment.name):
        session.add(experiment)
        session.flush()

def _upsert_gene(session, gene):

    if Gene.is_unique(session, gene.regionID, gene.sourceID, gene.name):
        session.add(gene)
        session.flush()

def _upsert_region(session, region):

    if Region.is_unique(session, region.chrom, region.start, region.end, region.strand):
        session.add(region)
        session.flush()

def _upsert_sample(session, sample):
    """
    @classmethod
    def is_unique(cls, session, name, treatment, cell_line, cancer):
    """
    if Sample.is_unique(session, sample.name, sample.treatment, sample.cell_line, sample.cancer):
        session.add(sample)
        session.flush()

def _upsert_source(session, source):

    if Source.is_unique(session, source.name):
        session.add(source)
        session.flush()

def _upsert_tf(session, tf):

    if TFBinding.is_unique(session, tf.regionID, tf.sampleID, tf.experimentID, tf.sourceID, tf.tf):
        session.add(tf)
        session.flush()

#--------------#
# Multiprocess #
#--------------#

def _process_data_in_chunks(data_file, insert_function, threads=1):

    from itertools import islice
    from multiprocessing import Pool

    # Initialize
    pool = Pool()

    # Get iterable
    iterable = GUDglobals.parse_tsv_file(data_file)

    # Get chunks
    chunks = _grouper(iterable)

    while True:

        # Groups chunks for multiprocessing
        grouped_chunks = [list(chunk) for chunk in islice(chunks, threads)]

        if grouped_chunks:
            pool.map(insert_function, grouped_chunks)

        else:
            break

        break

    # Close pool
    pool.close()

def _grouper(iterable, n=1000, fillvalue=None):

    import sys
    # Python 3+
    if sys.version_info > (3, 0):
        from itertools import zip_longest
    # Python 2.7
    else:
        from itertools import izip_longest as zip_longest

    args = [iter(iterable)] * n

    return(zip_longest(*args, fillvalue=fillvalue))