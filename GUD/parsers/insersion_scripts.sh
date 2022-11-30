###################### hg19 ######################
#clinvar
python -m GUD.parsers.clinvar2gud --genome hg19 \
--source_name ClinVar --clinvar_file \
/space/www/GUD_TABLES/Transfer_GRCh37_191007/clinvar.GRCh37.norm.anno.vcf \
-d hg19 -u gud_w

#strG
python -m GUD.parsers.str2gud --genome hg19 --source_name GangSTR --str_file \
/space/www/GUD_TABLES/Transfer_GRCh37_191007/Genomic_STR.GUDformatted.tsv \
--based 1 -d hg19 -u gud_w 

#strP
python -m GUD.parsers.str2gud --genome hg19 --source_name Richmond_github_20191007 --str_file \
/space/www/GUD_TABLES/Transfer_GRCh37_191007/Pathogenic_STR.GUDformatted.tsv \
--based 0 -d hg19 -u gud_w 

#CNV
python -m GUD.parsers.cnv2gud --genome hg19 --source_name dbVar --cnv_file \
/space/www/GUD_TABLES/Transfer_GRCh37_191007/GRCh37.nr_deletions.GUDformatted.tsv \
-d hg19 -u gud_w 

python -m GUD.parsers.cnv2gud --genome hg19 --source_name dbVar --cnv_file \
/space/www/GUD_TABLES/Transfer_GRCh37_191007/GRCh37.nr_duplications.GUDformatted.tsv \
-d hg19 -u gud_w 

###################### hg38 ######################
#clinvar
python -m GUD.parsers.clinvar2gud --genome hg38 \
--source_name ClinVar --clinvar_file \
/space/www/GUD_TABLES/Transfer_GRCh38_191007/clinvar.GRCh38.norm.anno.vcf \
-d hg38 -u gud_w 

#strG
python -m GUD.parsers.str2gud --genome hg38 --source_name GangSTR --str_file \
/space/www/GUD_TABLES/Transfer_GRCh38_191007/Genomic_STR.GUDformatted.tsv \
--based 1 -d hg38 -u gud_w 

#strP
python -m GUD.parsers.str2gud --genome hg38 --source_name Richmond_github_20191007 --str_file \
/space/www/GUD_TABLES/Transfer_GRCh38_191007/Pathogenic_STR.GUDformatted.tsv \
--based 0 -d hg38 -u gud_w 

#CNV
python -m GUD.parsers.cnv2gud --genome hg38 --source_name dbVar --cnv_file \
/space/www/GUD_TABLES/Transfer_GRCh38_191007/GRCh38.nr_deletions.GUDformatted.tsv \
-d hg38 -u gud_w 

python -m GUD.parsers.cnv2gud --genome hg38 --source_name dbVar --cnv_file \
/space/www/GUD_TABLES/Transfer_GRCh38_191007/GRCh38.nr_duplications.GUDformatted.tsv \
-d hg38 -u gud_w 