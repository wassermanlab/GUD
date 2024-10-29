##refgene/initialization
python -m GUD.parsers.refseq2gud --genome hg38 -d grch38 -u gud_w -p genebreaker --version "15-07-2020"

##clinvar (1based)
python -m GUD.parsers.clinvar2gud --genome hg38 \
--source_name ClinVar --clinvar_file \
/space/data_tables/Transfer_GRCh38_200326/clinvar.norm.anno.vcf \
-d grch38 -u gud_w -p genebreaker -P 3306 

##strG (1based)
nohup python -m GUD.parsers.str2gud --genome hg38 --source_name GangSTR --str_file \
/space/data_tables/Transfer_GRCh38_200326/Genomic_STR.GUDformatted.tsv \
--based 1 -d grch38 -u gud_w -p genebreaker -P 3306 &

##strP (0based)
python -m GUD.parsers.str2gud --genome hg38 --source_name Richmond_github_20191007 --str_file \
/space/data_tables/Transfer_GRCh38_200326/Pathogenic_STR.GUDformatted.tsv \
--based 0 -d grch38 -u gud_w -p genebreaker -P 3306

##CNV (0based)
python -m GUD.parsers.cnv2gud --genome hg38 --source_name dbVar --cnv_file \
/space/data_tables/Transfer_GRCh38_200326/GRCh38.nr_deletions.GUDformatted.tsv \
-d grch38 -u gud_w -p genebreaker -P 3306

##CNV (0based)
python -m GUD.parsers.cnv2gud --genome hg38 --source_name dbVar --cnv_file \
/space/data_tables/Transfer_GRCh38_200326/GRCh38.nr_duplications.GUDformatted.tsv \
-d grch38 -u gud_w -p genebreaker -P 3306

###################### hg19 ######################
##refgene/initialization

python -m GUD.parsers.refseq2gud --genome hg19 -d grch37 -u gud_w -p genebreaker --version "15-07-2020"

##clinvar (1based)
nohup python -m GUD.parsers.clinvar2gud --genome hg19 \
--source_name ClinVar --clinvar_file \
/space/data_tables/Transfer_GRCh37_200327/clinvar.norm.anno.vcf \
-d grch37 -u gud_w -p genebreaker -P 3306 &

##strG (1based)
nohup python -m GUD.parsers.str2gud --genome hg19 --source_name GangSTR --str_file \
/space/data_tables/Transfer_GRCh37_200327/Genomic_STR.GUDformatted.tsv \
--based 1 -d grch37 -u gud_w -p genebreaker -P 3306 &

##strP (0based)
nohup python -m GUD.parsers.str2gud --genome hg19 --source_name Richmond_github_20191007 --str_file \
/space/data_tables/Transfer_GRCh37_200327/Pathogenic_STR.GUDformatted.tsv \
--based 0 -d grch37 -u gud_w -p genebreaker -P 3306 &

##CNV (0based)
nohup python -m GUD.parsers.cnv2gud --genome hg19 --source_name dbVar --cnv_file \
/space/data_tables/Transfer_GRCh37_200327/GRCh37.nr_deletions.GUDformatted.tsv \
-d grch37 -u gud_w -p genebreaker -P 3306 &

##CNV (0based)
nohup python -m GUD.parsers.cnv2gud --genome hg19 --source_name dbVar --cnv_file \
/space/data_tables/Transfer_GRCh37_200327/GRCh37.nr_duplications.GUDformatted.tsv \
-d grch37 -u gud_w -p genebreaker -P 3306 &

