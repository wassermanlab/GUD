from GUD import GUDUtils
from GUD.api.api_helpers import set_db
from GUD.ORM import (Gene, ShortTandemRepeat, CNV, ClinVar, Conservation, CpGIsland,
                     DNAAccessibility, Enhancer, HistoneModification, RepeatMask, TAD, 
                     TFBinding, TSS, Chrom, Sample, Experiment, Source, Expression) 
import pandas as pd

switch = {
        "chroms": Chrom(),
        "clinvar": ClinVar(), 
        "copy_number_variants": CNV(),
        "conservation": Conservation(),
        "cpg_islands": CpGIsland(),
        "dna_accessibility": DNAAccessibility(),
        "enhancers": Enhancer(),
        "experiments": Experiment(),
        # "expression": Experession(),
        "genes": Gene(),
        "histone_modifications": HistoneModification(),
        "samples": Sample(),
        "short_tandem_repeats": ShortTandemRepeat(),
        "sources": Source(),
        "rmsk": RepeatMask(),
        "tads": TAD(),
        "tf_binding": TFBinding(),
        "tss": TSS()
    }

def is_gf1(table):
    """returns True if table is gf1 and False if not"""
    if table in ['clinvar', 'conservation', 'copy_number_variants', 'cpg_islands', 'genes', 'rmsk', 'short_tandem_repeats']:
        return True
    return False

def is_gf2(table):
    """returns True if table is gf2 and False if not"""
    if table in ['dna_accessibility','histone_modifications', 'tads', 'tf_binding']:
        return True
    return False

def get_unique_sources(table, session):
    if is_gf1(table): 
        resource = switch[table]
        return resource.get_unique_source_names(session)
    return None

def get_unique_samples(table, session):
    if is_gf2(table): 
        resource = switch[table]
        return resource.get_unique_sample_names(session)
    return None

def get_unique_experiments(table, session):
    if is_gf2(table): 
        resource = switch[table]
        return resource.get_unique_experiment_names(session)
    return None

def output_stats():

    set_db("hg38")
    engine_hg38, Session_hg38 = GUDUtils.get_engine_session(GUDUtils._get_db_name())
    session = Session_hg38()

    tables = []
    mb = []
    rows = []
    unique_source_names = []
    unique_sample_names = []
    unique_experiment_names = []


    with engine_hg38.connect() as con:
        rs = con.execute('''
            SELECT table_name, (data_length+index_length)/power(1024,2) tablesize_mb, table_rows
            FROM information_schema.tables
            WHERE table_schema='hg38';
        ''')

    for row in rs: 
        tables.append(row[0])
        mb.append(row[1])
        rows.append(row[2])

    for t in tables:
        unique_source_names.append(get_unique_sources(t, session))
        unique_sample_names.append(get_unique_samples(t, session))
        unique_experiment_names.append(get_unique_experiments(t, session)) 

    # output experiments table
    # output samples table
    # output sources table

    session.close()
    engine_hg38.dispose() 

output_stats()
# # experiment 
# total_rows = session.query(Experiment).count()
# distinct_names_count = session.query(Experiment.name).distinct(Experiment.name).group_by(Experiment.name).count()
# d = {"total_rows": [total_rows], "distinct_names_count": [distinct_names_count]}
# df = pd.DataFrame(data=d)
# df.to_csv('GUD/api/static/stats/experiments.csv')
# all_rows = session.query(Experiment).all()
# header = Experiment.__table__.columns.keys()
# df = pd.DataFrame(columns=header)
# for i in range(len(all_rows)):
#    row = all_rows[i]
#    df.loc[i] = [row.uid, row.name, row.experiment_metadata, row.metadata_descriptor]
# df.to_csv('GUD/api/static/stats/all_experiments.csv')

# # sample 
# total_rows = session.query(Sample).count()
# distinct_names_count = session.query(Sample.name).distinct(Sample.name).group_by(Sample.name).count()
# d = {"total_rows": [total_rows], "distinct_names_count": [distinct_names_count]}
# df = pd.DataFrame(data=d)
# df.to_csv('GUD/api/static/stats/samples.csv')
# all_rows = session.query(Sample).all()
# header = Sample.__table__.columns.keys()
# df = pd.DataFrame(columns=header)
# for i in range(len(all_rows)):
#    row = all_rows[i]
#    df.loc[i] = [row.uid, row.name, row.X, row.Y, row.treatment, row.cell_line, row.cancer]
# df.to_csv('GUD/api/static/stats/all_samples.csv')

# # sources 
# total_rows = session.query(Source).count()
# distinct_names_count = session.query(Source.name).distinct(Source.name).group_by(Source.name).count()
# d = {"total_rows": [total_rows], "distinct_names_count": [distinct_names_count]}
# df = pd.DataFrame(data=d)
# df.to_csv('GUD/api/static/stats/source.csv')
# all_rows = session.query(Source).all()
# header = Source.__table__.columns.keys()
# df = pd.DataFrame(columns=header)
# for i in range(len(all_rows)):
#    row = all_rows[i]
#    df.loc[i] = [row.uid, row.name, row.source_metadata, row.metadata_descriptor, row.url, row.insert_date]
# df.to_csv('GUD/api/static/stats/all_sources.csv')

# #chrom
# total_rows = session.query(Chrom).count()
# d = {"total_rows": [total_rows]}
# df = pd.DataFrame(data=d)
# df.to_csv('GUD/api/static/stats/chrom.csv')

# #clinvar
# total_rows = session.query(ClinVar).count()
# d = {"total_rows": [total_rows]}
# df = pd.DataFrame(data=d)
# df.to_csv('GUD/api/static/stats/clinvar.csv')

# #conservation
# #copy_number_variant
# #cpg_island
# #dna_accessibility
# #enhancer
# #expression
# #gene
# #histone_modification
# #repeat_mask
# #short_tandem_repeat
# #tad
# #tf_binding 
# #tss



