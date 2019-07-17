from GUD.api import app
from flask import request, jsonify, render_template, url_for
from werkzeug.exceptions import HTTPException, NotFound, BadRequest
from GUD.api.routes_api import * 
import json
import sys, os

@app.route('/')
def home():
    return render_template('home.html')

@app.route('/live_api')
def live_api():
    filename = os.path.join(app.static_folder, 'docs.json')
    data = None
    with open(filename) as json_file:
        data = json.load(json_file)
    routes = data.keys()
    return render_template('live_api.html', routes = routes)

@app.route('/docs')
def docs():
    return "DOCS"
    # return render_template('home.html')

@app.route('/contact')
def contact():
    return "CONTACT"
    # return render_template('home.html')