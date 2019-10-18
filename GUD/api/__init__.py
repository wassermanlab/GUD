from flask import Flask
from flask_limiter import Limiter
from flask_limiter.util import get_remote_address
from GUD import GUDUtils


app = Flask(__name__)
app.config.from_mapping(
    SECRET_KEY='dev',
)
limiter = Limiter(
    app,
    key_func=get_remote_address,
    default_limits=["5 per second"]
)
import GUD.api.routes

