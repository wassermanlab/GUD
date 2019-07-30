###################### hg19 ######################
#clinvar
python -m GUD.parsers.clinvar2gud --genome hg19 \
--source_name ClinVar --clinvar_file \
/space/home/tavshalom/GUD_TABLES/Transfer_GRCh37_190723/clinvar.GRCh37.norm.anno.vcf \
-d hg19 -u gud_w --threads 6

#strG
python -m GUD.parsers.str2gud --genome hg19 --source_name GangSTR --str_file \
/space/home/tavshalom/GUD_TABLES/Transfer_GRCh37_190723/Genomic_STR.GUDformatted.tsv \
--based 1 -d hg19 -u gud_w --threads 6

#strP
python -m GUD.parsers.str2gud --genome hg19 --source_name Richmond_github_20190723 --str_file \
/space/home/tavshalom/GUD_TABLES/Transfer_GRCh37_190723/Pathogenic_STR.GUDformatted.tsv \
--based 0 -d hg19 -u gud_w --threads 6

#CNV
python -m GUD.parsers.cnv2gud --genome hg19 --source_name dbVar --cnv_file \
/space/home/tavshalom/GUD_TABLES/Transfer_GRCh37_190723/GRCh37.nr_deletions.GUDformatted.tsv \
-d hg19 -u gud_w --threads 6

python -m GUD.parsers.cnv2gud --genome hg19 --source_name dbVar --cnv_file \
/space/home/tavshalom/GUD_TABLES/Transfer_GRCh37_190723/GRCh37.nr_duplications.GUDformatted.tsv \
-d hg19 -u gud_w --threads 6

###################### hg38 ######################
#clinvar
python -m GUD.parsers.clinvar2gud --genome hg38 \
--source_name ClinVar --clinvar_file \
/space/home/tavshalom/GUD_TABLES/Transfer_GRCh38_190723/clinvar.GRCh38.norm.anno.vcf \
-d hg38 -u gud_w --threads 6

#strG
python -m GUD.parsers.str2gud --genome hg38 --source_name GangSTR --str_file \
/space/home/tavshalom/GUD_TABLES/Transfer_GRCh38_190723/Genomic_STR.GUDformatted.tsv \
--based 1 -d hg38 -u gud_w --threads 6

#strP
python -m GUD.parsers.str2gud --genome hg38 --source_name Richmond_github_20190723 --str_file \
/space/home/tavshalom/GUD_TABLES/Transfer_GRCh38_190723/Pathogenic_STR.GUDformatted.tsv \
--based 0 -d hg38 -u gud_w --threads 6

#CNV
python -m GUD.parsers.cnv2gud --genome hg38 --source_name dbVar --cnv_file \
/space/home/tavshalom/GUD_TABLES/Transfer_GRCh38_190723/GRCh38.nr_deletions.GUDformatted.tsv \
-d hg38 -u gud_w --threads 6 

python -m GUD.parsers.cnv2gud --genome hg38 --source_name dbVar --cnv_file \
/space/home/tavshalom/GUD_TABLES/Transfer_GRCh38_190723/GRCh38.nr_duplications.GUDformatted.tsv \
-d hg38 -u gud_w --threads 6 