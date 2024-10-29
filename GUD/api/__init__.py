from flask import Flask
from GUD import GUDUtils
from GUD.api.api_helpers import set_db
from werkzeug.exceptions import BadRequest

app = Flask(__name__)
app.config.from_pyfile('config.py')
set_db("grch37")
engine_grch37, Session_grch37 = GUDUtils.get_engine_session(GUDUtils._get_db_name())

set_db("grch38")
engine_grch38, Session_grch38 = GUDUtils.get_engine_session(GUDUtils._get_db_name())

def get_engine_session(db): 
    if db == "grch37":
        return engine_grch37, Session_grch37()
    elif db == "grch38":
        return engine_grch38, Session_grch38()
    else:
        raise BadRequest('database must be grch37 or grch38.')

import GUD.api.routes

