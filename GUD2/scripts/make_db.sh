#DB initialization (chrom, region, sample, source, experiment tables)
/usr/local/python2.7.14/bin/python -m GUD2.scripts.initialize_gud_db hg19 -u ontarget_w -d tamar_test

#str insertion 
nohup /usr/local/python2.7.14/bin/python -m GUD2.scripts.str2gud /space/data/gangSTR/hg19_ver10.sorted.bed -u ontarget_w -d tamar_test &

# insert the gene table
nohup /usr/local/python2.7.14/bin/python -m GUD2.scripts.refgene2gud hg19 -u ontarget_w -d tamar_test &

#conservation
nohup /usr/local/python2.7.14/bin/python -m GUD2.scripts.conservation2gud hg19 -u ontarget_w -d tamar_test &

#rmsk
nohup /usr/local/python2.7.14/bin/python -m GUD2.scripts.rmsk2gud hg19 -u ontarget_w -d tamar_test &

#fantom2gud
nohup /usr/local/python2.7.14/bin/python -m GUD2.scripts.fantom2gud hg19 -u ontarget_w -d tamar_test &
