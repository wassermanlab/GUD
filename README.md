# Genomic Unification Database
test
> This is a database which centralizes and unifies genomic data in a universal manner for specific reference datasets.

## Manifest
+ GUD/ORM - Object relational mapping classes
+ GUD/parsers - Parsers for upserting data into GUD
+ GUD/scripts - General scripts

## Requirements
GUD requires the following dependencies:
* [`GNU Core Utilities`](https://www.gnu.org/software/coreutils/)
* [`MySQL`](https://www.mysql.com)
* [`Python`](https://www.python.org) `≥2.7` or `3.x` (recommended) with:
    - [`Biopython`](https://biopython.org)
    - [`Core API`](https://www.coreapi.org/)
    - [`flask`](https://flask.palletsprojects.com/en/1.0.x/)
    - [`Flask-SQLAlchemy`](https://flask-sqlalchemy.palletsprojects.com/en/2.x/)
    - [`FuzzyWuzzy`](https://chairnerd.seatgeek.com/fuzzywuzzy-fuzzy-string-matching-in-python/)
    - [`interval-binning`](https://interval-binning.readthedocs.io/en/latest/)
    - [`NumPy`](https://numpy.org/)
    - [`pandas`](https://pandas.pydata.org/)
    - [`pybedtools`](https://daler.github.io/pybedtools/)
    - [`PyMySQL`](https://pymysql.readthedocs.io/en/latest/)
    - [`pyliftover`](https://github.com/konstantint/pyliftover)
    - [`requests`](https://requests.readthedocs.io/en/master/)
    - [`SQLAlchemy`](https://www.sqlalchemy.org)
    - [`SQLAlchemy-FullText-Search`](https://github.com/mengzhuo/sqlalchemy-fulltext-search)
    - [`SQLAlchemy-Utils`](https://sqlalchemy-utils.readthedocs.io/en/latest/)

## Configuration

```
conda create -n gud -c bioconda -c conda-forge python=3.7 biopython coreutils flask \
    flask-sqlalchemy fuzzywuzzy pandas pymysql pyliftover python-coreapi requests
pip install flask-limiter interval-binning python-Levenshtein SQLAlchemy-FullText-Search \
    sqlalchemy-utils warnings
```

## START UP SERVER

```
export FLASK_APP=GUD/api
export FLASK_ENV=development
flask run
```
