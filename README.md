# Genomic Unification Database
test
> This is a database which centralizes and unifies genomic data in a universal manner for specific reference datasets.

## Manifest
+ GUD/ORM - Object relational mapping classes
+ GUD/parsers - Parsers for upserting data into GUD
+ GUD/scripts - General scripts

## Requirements
GUD requires the following dependencies:
* [`MySQL`](https://www.mysql.com)
* [`Parallel`](https://www.gnu.org/software/parallel/)
* [`Python`](https://www.python.org) `≥2.7` or `3.x` with:
    - [`Biopython`](https://biopython.org)
    - [`interval-binning`](https://interval-binning.readthedocs.io/en/latest/)
    - [`macs2`](https://github.com/taoliu/MACS/)
    - [`pandas`](https://pandas.pydata.org/)
    - [`PyMySQL`](https://pymysql.readthedocs.io/en/latest/)
    - [`SQLAlchemy`](https://www.sqlalchemy.org)
    - `SQLAlchemy-FullText-Search`
    - `sqlalchemy-utils`

## INSTALLATION

```bash
conda create -n gud -c bioconda python=3.7 pybedtools biopython pymysql
pip install interval-binning SQLAlchemy-FullText-Search sqlalchemy-utils
```

## SETUP

```bash
conda env create -f environment.yml
```

## START UP SERVER OLD

```bash
conda activate GUD
export FLASK_APP=GUD/api
export FLASK_ENV=development
```

## START UP SERVER

```
FLASK_APP=GUD/api FLASK_ENV=development python -m flask run

FLASK_APP=GUD/api FLASK_ENV=development flask run

export FLASK_APP=GUD/api
export FLASK_ENV=development
flask run
```
