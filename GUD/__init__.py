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

from sqlalchemy import create_engine
from sqlalchemy.orm import scoped_session, sessionmaker

__all__ = ["ORM", "parsers"]

class GUDUtilities:
    """
    Contains functions designed to work through the entire module.
    """

    #--------------#
    # Defaults     #
    #--------------#

    def __init__(self):
        """
        MySQL options:
        @param user = User for login (i.e. option "-u")
        @param pwd  = Password to use when connecting (i.e. option "-p")
        @param host = Server host (i.e. option "-h")
        @param port = Port number to use for connection (i.e. option "-P")
        @param db   = Database to use (i.e. option "-D")
        """

        # i.e. defaults for GUD at the Wasserman Lab
        self._user = "gud_r"
        self._pwd  = ""
        self._host = "gud.cmmt.ubc.ca"
        self._port = 5506
        self._db   = "hg38"

    @property
    def user(self):
        """
        MySQL user:
        @rtype = {String}
        """
        return(self._user)

    @user.setter
    def user(self, value):
        self._user = str(value)

    @property
    def pwd(self):
        """
        MySQL password:
        @rtype = {String}
        """
        return(self._pwd)

    @pwd.setter
    def pwd(self, value):
        self._pwd = str(value)

    @property
    def host(self):
        """
        MySQL host:
        @rtype = {String}
        """
        return(self._host)

    @host.setter
    def host(self, value):
        self._host = str(value)

    @property
    def port(self):
        """
        MySQL port (e.g. "3306"):
        @rtype = {Integer}
        """
        return(self._port)

    @port.setter
    def port(self, value):
        self._port = int(value)

    @property
    def db(self):
        """
        MySQL database:
        @rtype = {String}
        """
        return(self._db)

    @db.setter
    def db(self, value):
        self._db = str(value)

    #--------------#
    # SQLalchemy   #
    #--------------#

    def get_engine(self):
        """
        Create an SQLAlchemy {Engine} and bind it to a {Session} factory:
        @rtype = {Session}
        """

        db_name = self._get_db_name()

        try:
            engine, Session = self._get_engine_session(db_name)
        except:
            raise ValueError("Could not connect to GUD: %s" % db_name)

        return(engine)

    def get_session(self):
        """
        Create an SQLAlchemy {Engine} and bind it to a {Session} factory:
        @rtype = {Session}
        """

        db_name = self._get_db_name()

        try:
            engine, Session = self._get_engine_session(db_name)
        except:
            raise ValueError("Could not connect to GUD: %s" % db_name)

        return(Session)

    def _get_db_name(self):
        return("mysql+pymysql://{}:{}@{}:{}/{}".format(self.user, self.pwd, self.host, self.port, self.db))

    def _get_engine_session(self, db_name):

        # Initialize
        engine = create_engine(db_name, pool_pre_ping=True, pool_size=20, max_overflow=0)
        session_factory = sessionmaker(bind=engine)
        Session = scoped_session(session_factory)

        return(engine, Session)

GUDUtils = GUDUtilities()