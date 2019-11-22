from flask import Flask
from flask_limiter import Limiter
from flask_limiter.util import get_remote_address
from GUD import GUDUtils
from GUD.api.api_helpers import set_db
from werkzeug.exceptions import BadRequest

app = Flask(__name__)
# app.config.from_mapping(SECRET_KEY='dev',)
app.config.from_pyfile('config.py')
limiter = Limiter(
    app,
    key_func=get_remote_address,
    default_limits=["5 per second"]
)
set_db("hg19")
engine_hg19, Session_hg19 = GUDUtils.get_engine_session(GUDUtils._get_db_name())

set_db("hg38")
engine_hg38, Session_hg38 = GUDUtils.get_engine_session(GUDUtils._get_db_name())

set_db("test")
engine_test, Session_test = GUDUtils.get_engine_session(GUDUtils._get_db_name())

set_db("test_hg38_chr22")
engine_test_hg38_chr22, Session_test_hg38_chr22 = GUDUtils.get_engine_session(GUDUtils._get_db_name())

def get_engine_session(db): 
    if db == "hg19":
        return engine_hg19, Session_hg19()
    elif db == "hg38":
        return engine_hg38, Session_hg38()
    elif db == "test":
        return engine_test, Session_test()
    elif db == "test_hg38_chr22":
        return engine_test_hg38_chr22, Session_test_hg38_chr22()
    else:
        raise BadRequest('database must be hg19 or hg38 or test or test_hg38_chr22')

import GUD.api.routes

