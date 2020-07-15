
import csv
from sqlalchemy import create_engine, update
from sqlalchemy.orm import (
    scoped_session,
    sessionmaker
)
from sqlalchemy_utils import database_exists

# Import from GUD module
from GUD.ORM.sample import Sample
from GUD import GUDUtils
from GUD.parsers import ParseUtils

# Get database name
db_name = GUDUtils._get_db_name()

# Get engine/session
engine, Session = GUDUtils.get_engine_session(db_name)

# Initialize parser utilities
ParseUtils.genome = "hg38"
ParseUtils.dbname = db_name
ParseUtils.engine = engine

# Start a new session
session = Session()

with open('samples.csv', 'r') as file:
    reader = csv.reader(file)
    row_count = 0
    for row in reader:
        if row_count == 0:
            row_count += 1
        else:
            row_count += 1
            
            session.query(Sample).filter(Sample.uid == int(row[0])).\
            update({
                Sample.name: row[1],
                Sample.X: row[2],
                Sample.Y: row[3],
                Sample.treatment: row[4],
                Sample.cell_line: row[5],
                Sample.cancer: row[6],})
            session.commit()

# This is ABSOLUTELY necessary to prevent MySQL from crashing!
session.close()
engine.dispose()