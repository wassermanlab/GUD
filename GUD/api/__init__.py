from flask import Flask
from flask_sqlalchemy import SQLAlchemy

app = Flask(__name__)
app.config.from_mapping(
    SECRET_KEY='dev',
)
# app.config['SQLALCHEMY_DATABASE_URI'] ="mysql://ontarget_r:@ontarget.cmmt.ubc.ca:5506/tamar_test"
# app.config['SQLALCHEMY_DATABASE_URI'] ="mysql://ontarget_r:@ontarget.cmmt.ubc.ca:5506/tamar_test"
# db = SQLAlchemy(app)

import GUD.api.routes
