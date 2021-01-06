coverage run -m -a GUD.tests.test_chrom
coverage run -m -a GUD.tests.test_conservation
coverage run -m -a GUD.tests.test_gene
coverage run -m -a GUD.tests.test_api_chrom
coverage run -m -a GUD.tests.test_api_conservation
coverage run -m -a GUD.tests.test_api_gene
coverage report GUD/ORM/*.py
coverage report GUD/api/*.py
coverage html