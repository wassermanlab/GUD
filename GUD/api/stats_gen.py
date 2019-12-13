from GUD import GUDUtils
from GUD.api.api_helpers import set_db
from GUD.ORM import (Gene, ShortTandemRepeat, CNV, ClinVar, Conservation, CpGIsland,
                     DNAAccessibility, Enhancer, HistoneModification, RepeatMask, TAD,
                     TFBinding, TSS, Chrom, Sample, Experiment, Source, Expression)
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm

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
    if table in ['clinvar', 'conservation', 'copy_number_variants', 'cpg_islands', 'genes', 'rmsk', 'short_tandem_repeats', 'dna_accessibility', 'histone_modifications', 'tads', 'tf_binding']:
        return True
    return False


def is_gf2(table):
    """returns True if table is gf2 and False if not"""
    if table in ['dna_accessibility', 'histone_modifications', 'tads', 'tf_binding']:
        return True
    return False


def get_unique_sources(table, session):
    if is_gf1(table):
        resource = switch[table]
        return ";".join(resource.get_unique_source_names(session))
    return None


def get_unique_samples(table, session):
    if is_gf2(table):
        resource = switch[table]
        return ";".join(resource.get_unique_sample_names(session))
    return None


def get_unique_experiments(table, session):
    if is_gf2(table):
        resource = switch[table]
        return ";".join(resource.get_unique_experiment_names(session))
    return None


def output_stats():

    set_db("hg38")
    engine_hg38, Session_hg38 = GUDUtils.get_engine_session(
        GUDUtils._get_db_name())
    session = Session_hg38()

    tables = []
    mb = []
    percent = []
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

    total = sum(mb)
    percent = [(i/total)*100 for i in mb]

    # output summary
    d = {"table": tables, "mb": mb, "percent": percent,
         "rows": rows, "unique_source_names": unique_source_names,
         "unique_sample_names": unique_sample_names, "unique_experiment_names": unique_experiment_names}
    df = pd.DataFrame(data=d)
    df.to_csv('GUD/api/static/stats/hg38_summary_stats.csv', index=False)

    # output pie chart 
    labels = tables
    sizes = percent
    theme = plt.get_cmap('tab20c')
    colors = ([theme(1. * i / len(sizes))for i in range(len(sizes))])
    patches, texts = plt.pie(sizes, colors=colors,)
    plt.legend(patches, labels=['%s, %1.1f %%' % (l, s) for l, s in zip(labels, sizes)], loc="center left",bbox_to_anchor=(-0.30, .5))
    plt.savefig('GUD/api/static/stats/hg38_piechart.png')

    # output experiments table
    all_rows = session.query(Experiment).all()
    header = Experiment.__table__.columns.keys()
    df = pd.DataFrame(columns=header)
    for i in range(len(all_rows)):
        row = all_rows[i]
        df.loc[i] = [row.uid, row.name,
                     row.experiment_metadata, row.metadata_descriptor]
    df.to_csv('GUD/api/static/stats/hg38_all_experiments.csv', index=False)

    # output samples table
    all_rows = session.query(Sample).all()
    header = Sample.__table__.columns.keys()
    df = pd.DataFrame(columns=header)
    for i in range(len(all_rows)):
        row = all_rows[i]
        df.loc[i] = [row.uid, row.name, row.X, row.Y,
                     row.treatment, row.cell_line, row.cancer]
    df.to_csv('GUD/api/static/stats/hg38_all_samples.csv', index=False)

    # output sources table
    all_rows = session.query(Source).all()
    header = Source.__table__.columns.keys()
    df = pd.DataFrame(columns=header)
    for i in range(len(all_rows)):
        row = all_rows[i]
        df.loc[i] = [row.uid, row.name, row.source_metadata,
                     row.metadata_descriptor, row.url, row.insert_date]
    df.to_csv('GUD/api/static/stats/hg38_all_sources.csv', index=False)

    session.close()
    engine_hg38.dispose()


output_stats()


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
