##DB initialization (chrom, region, sample, source, experiment tables)
/usr/local/python2.7.14/bin/python -m GUD.parsers.initialize --genome hg19 -u ontarget_w -d tamar_test

##str insertion 
nohup /usr/local/python2.7.14/bin/python -m GUD.parsers.str2gud /space/data/gangSTR/hg19_ver10.sorted.0based.bed -u ontarget_w -d tamar_test &

##pathogenic str insertion 
/usr/local/python2.7.14/bin/python -m GUD.parsers.path_str2gud ../Pathogenic_STR_Table_Shifted.bed -u ontarget_w -d tamar_test --source PathogenicSTR-20181219

## insert the gene table
nohup /usr/local/python2.7.14/bin/python -m GUD.parsers.refgene2gud --genome hg19 -u ontarget_w -d tamar_test &

##clinvar insertion 
nohup /usr/local/python2.7.14/bin/python -m GUD.parsers.clinvar2gud ../clinvar.norm.anno.vcf -u ontarget_w -d tamar_test &

##cnv insertion 
nohup /usr/local/python2.7.14/bin/python -m GUD.parsers.cnv2gud ../isca_Table.tsv -u ontarget_w -d tamar_test &

# #conservation
# nohup /usr/local/python2.7.14/bin/python -m GUD2.scripts.conservation2gud hg19 -u ontarget_w -d tamar_test &

# #rmsk
# nohup /usr/local/python2.7.14/bin/python -m GUD2.scripts.rmsk2gud hg19 -u ontarget_w -d tamar_test &

# #fantom2gud
# nohup /usr/local/python2.7.14/bin/python -m GUD2.scripts.fantom2gud hg19 -u ontarget_w -d tamar_test &

# #clinvar insertion 
# nohup /usr/local/python2.7.14/bin/python -m GUD2.scripts.clinvar2gud ../clinvar.norm.anno.vcf -u ontarget_w -d tamar_test &

# #cnv insertion 
# nohup /usr/local/python2.7.14/bin/python -m GUD2.scripts.cnv2gud <file> -u ontarget_w -d tamar_test &
