#DB initialization (chrom, region, sample, source, experiment tables)
/usr/local/python2.7.14/bin/python -m GUD2.scripts.initialize_gud_db hg19 -u ontarget_w -d $DB

#str insertion 
/usr/local/python2.7.14/bin/python -m GUD2.scripts.str2gud /space/data/gangSTR/hg19_ver10.sorted.bed -u ontarget_w -d $DB