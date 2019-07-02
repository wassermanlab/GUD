# Genomic Universal Database

> This is a database which centralizes and unifies genomic data in a universal manner for specific reference datasets.

## Manifest
+ GUD/ORM - Object relational mapping classes
+ GUD/parsers - Parsers for upserting data into GUD
+ GUD/scripts - General scripts

## Requirements
GUD requires the following dependencies:
* [`MySQL`](https://www.mysql.com)
* [`Python`](https://www.python.org) `≥2.7` or `3.x` with:
    - [`Biopython`](https://biopython.org)
    - [`interval-binning`](https://interval-binning.readthedocs.io/en/latest/)
    - [`pandas`](https://pandas.pydata.org/)
    - [`PyMySQL`](https://pymysql.readthedocs.io/en/latest/)
    - [`SQLAlchemy`](https://www.sqlalchemy.org)
    - `SQLAlchemy-FullText-Search`
    - `sqlalchemy-utils`

```
conda create -n gud -c bioconda python=3.7 pybedtools biopython pymysql
pip install interval-binning SQLAlchemy-FullText-Search sqlalchemy-utils
```