coverage run -m -a GUD.tests.test_chrom
coverage run -m -a GUD.tests.test_clinvar
coverage run -m -a GUD.tests.test_cnv
coverage run -m -a GUD.tests.test_conservation
coverage run -m -a GUD.tests.test_gene
coverage run -m -a GUD.tests.test_str
coverage run -m -a GUD.tests.test_api_chrom
coverage run -m -a GUD.tests.test_api_clinvar
coverage run -m -a GUD.tests.test_api_cnv
coverage run -m -a GUD.tests.test_api_conservation
coverage run -m -a GUD.tests.test_api_gene
coverage run -m -a GUD.tests.test_api_str
coverage report GUD/ORM/*.py
coverage report GUD/api/*.py
coverage html