from GUD.api import app
from flask import request, jsonify
from werkzeug.exceptions import HTTPException, NotFound, BadRequest
from GUD.api.routes_api import * 

@app.route('/')
def index():
    return 'HOME'