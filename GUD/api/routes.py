from GUD.api import app
from flask import request, jsonify, render_template
from werkzeug.exceptions import HTTPException, NotFound, BadRequest
from GUD.api.routes_api import * 

@app.route('/')
def home():
    return render_template('home.html')

@app.route('/live_api')
def live_api():
    return render_template('live_api.html')

@app.route('/docs')
def docs():
    return "DOCS"
    # return render_template('home.html')

@app.route('/contact')
def contact():
    return "CONTACT"
    # return render_template('home.html')