from GUD.api import app
from flask import request, jsonify, render_template, url_for
from werkzeug.exceptions import HTTPException, NotFound, BadRequest
from GUD.api.routes_api import * 
import json
import sys, os, html

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
    with open(filename) as json_file:
        data = json.load(json_file)
    return render_template('live_api.html', routes = routes, resources = data, address_base = app.config["ADDRESS_BASE"])

@app.route('/docs')
def docs():
    filename = os.path.join(app.static_folder, 'docs.json')
    docs = None
    with open(filename) as json_file:
        docs = json.load(json_file)
    filename2 = os.path.join(app.static_folder, 'docs.titles.json')
    titles = None
    with open(filename2) as json_file:
        titles = json.load(json_file)
    return render_template('docs.html', docs = docs, titles = titles)

@app.route('/contact')
def contact():
    return render_template('contact.html')