#!/usr/bin/env python

import os, sys
from binning import containing_bins, contained_bins 
from sqlalchemy import create_engine, Index
from sqlalchemy.orm import Session
import warnings

# Append GUD module to path
module_path = os.path.join(os.path.dirname(__file__), os.pardir)
sys.path.append(os.path.join(module_path, os.pardir))

# Import from GUD/OnTarget
#from GUD.ORM.chrom_size import ChromSize
#from GUD.ORM.conservation import Conservation
#from GUD.ORM.dna_accessibility import DnaAccessibility
#from GUD.ORM.gene import Gene
#from GUD.ORM.histone_modification import HistoneModification
#from GUD.ORM.repeat_mask import RepeatMask
#from GUD.ORM.tad import Tad
#from GUD.ORM.tf_binding import TfBinding
from GUD.ORM.tss import TSS

# Establish a SQLalchemy session w/ GUD
db_name = "mysql://{}:@{}:{}/{}".format("ontarget_r",
    "ontarget.cmmt.ubc.ca", "5506", "hg19")
try:
    engine = create_engine(db_name, echo=False)
    session = Session(engine)
except:
    raise ValueError("Cannot connect to GUD: %s" % "hg19")

# Initialize
chrom = "chr10"
start = 100106388
end = 100211911
#histone_types = "H3K4me1,H3K4me2,H3K4me3,H3K9ac,H3K9me1,H3K9me3,H3K27ac,H3K27me3".split(",")
#repeat_classes = "LINE,SINE,LTR".split(",")
gene = "PYROXD2"
tss = 2
samples = [
    "endothelial cells (aorta)",
    "endothelial cells (artery)",
    "endothelial cells (lymph node)",
    "endothelial cells (microvasculature)",
    "endothelial cells (thorax)",
    "endothelial cells (umbilical vein)",
    "endothelial cells (vein)",
    "endothelial cells (liver sinusoid)",
    "endothelial cells (renal glomerulus)"
]

# Get all bins overlapping range
bins = set(containing_bins(start, end) + contained_bins(start, end))
print(bins)

## Get chromosome sizes
#chrom_sizes = ChromSize.chrom_sizes(session)
#print(chrom_sizes)
#
## Select stuff by bins
#conservation = Conservation.select_by_bin_range(session,
#    chrom, start, end, list(bins))
#print("\nConserved Regions:")
#for c in conservation: print(c)
#
#dna_accessibility = DnaAccessibility.select_by_bin_range(
#    session, chrom, start, end, bins=list(bins))
#print("\nDNA Accessibility:")
#for da in dna_accessibility: print(da)
#
#genes = Gene.select_by_bin_range(session, chrom, start,
#    end, bins=list(bins))
#print("\nGenes:")
#for g in genes: print(g)
#
#histone_modification = HistoneModification.select_by_bin_range(
#    session, chrom, start, end, histone_types=histone_types,
#    bins=list(bins))
#print("\nHistone Modifications:")
#for hm in histone_modification: print(hm)
#
#repeats = RepeatMask.select_by_bin_range(session, chrom,
#    start, end, repeat_classes=repeat_classes, bins=list(bins))
#print("\nRepeats:")
#for r in repeats: print(r)
#
#tads = Tad.select_by_bin_range(session, chrom, start,
#    end, bins=list(bins))
#print("\nTADs:")
#for tad in tads: print(tad)
#
#tf_binding = TfBinding.select_by_bin_range(session, chrom,
#    start, end, bins=list(bins))
#print("\nTF-Binding:")
#for tf in tf_binding: print(tf)

feats = TSS.select_by_bin_range(session, chrom, start, end,
    bins=list(bins))
print("\nTSS:")
for t in feats: print(t)

feats = TSS.select_by_gene(session, gene)
print("\n%s:" % gene)
for t in feats: print(t)

feats = TSS.select_by_gene_tss(session, gene, tss)
print("\n%s (%s):" % (gene, tss))
for t in feats: print(t)

feats = TSS.select_by_sample(session, sample=samples, min_tpm=10.0)
print("\nTSSs in endothelial cells:")
for t in feats: print(t)
