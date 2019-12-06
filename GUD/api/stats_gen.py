from GUD import GUDUtils
from GUD.api.api_helpers import set_db
from GUD.ORM import (Gene, ShortTandemRepeat, CNV, ClinVar, Conservation, CpGIsland,
                     DNAAccessibility, Enhancer, HistoneModification, RepeatMask, TAD, 
                     TFBinding, TSS, Chrom, Sample, Experiment, Source, Expression) 
import pandas as pd

set_db("hg38")
engine_hg38, Session_hg38 = GUDUtils.get_engine_session(GUDUtils._get_db_name())
session = Session_hg38()

# experiment 
total_rows = session.query(Experiment).count()
distinct_names_count = session.query(Experiment.name).distinct(Experiment.name).group_by(Experiment.name).count()
d = {"total_rows": [total_rows], "distinct_names_count": [distinct_names_count]}
df = pd.DataFrame(data=d)
df.to_csv('GUD/api/static/stats/experiments.csv')
all_rows = session.query(Experiment).all()
header = Experiment.__table__.columns.keys()
df = pd.DataFrame(columns=header)
for i in range(len(all_rows)):
   row = all_rows[i]
   df.loc[i] = [row.uid, row.name, row.experiment_metadata, row.metadata_descriptor]
df.to_csv('GUD/api/static/stats/all_experiments.csv')

# sample 
total_rows = session.query(Sample).count()
distinct_names_count = session.query(Sample.name).distinct(Sample.name).group_by(Sample.name).count()
d = {"total_rows": [total_rows], "distinct_names_count": [distinct_names_count]}
df = pd.DataFrame(data=d)
df.to_csv('GUD/api/static/stats/samples.csv')
all_rows = session.query(Sample).all()
header = Sample.__table__.columns.keys()
df = pd.DataFrame(columns=header)
for i in range(len(all_rows)):
   row = all_rows[i]
   df.loc[i] = [row.uid, row.name, row.X, row.Y, row.treatment, row.cell_line, row.cancer]
df.to_csv('GUD/api/static/stats/all_samples.csv')

# sources 
total_rows = session.query(Source).count()
distinct_names_count = session.query(Source.name).distinct(Source.name).group_by(Source.name).count()
d = {"total_rows": [total_rows], "distinct_names_count": [distinct_names_count]}
df = pd.DataFrame(data=d)
df.to_csv('GUD/api/static/stats/source.csv')
all_rows = session.query(Source).all()
header = Source.__table__.columns.keys()
df = pd.DataFrame(columns=header)
for i in range(len(all_rows)):
   row = all_rows[i]
   df.loc[i] = [row.uid, row.name, row.source_metadata, row.metadata_descriptor, row.url, row.insert_date]
df.to_csv('GUD/api/static/stats/all_sources.csv')

#chrom
#clinvar
#conservation
#copy_number_variant
#cpg_island
#dna_accessibility
#enhancer
#expression
#gene
#histone_modification
#repeat_mask
#short_tandem_repeat
#tad
#tf_binding 
#tss

session.close()
engine_hg38.dispose()  

