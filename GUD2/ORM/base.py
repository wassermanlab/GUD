from sqlalchemy import create_engine, ForeignKey
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column, Integer, String
from sqlalchemy import Table, Text

engine = create_engine('mysql://test:test@localhost/test1',
                    echo=False)

Base = declarative_base()