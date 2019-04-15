#!/usr/bin/env python

import hashlib
import os
import re

# Import from GUD
from GUD import GUDglobals
from GUD.ORM.enhancer import Enhancer
from GUD.ORM.sample import Sample
from GUD.ORM.tad import TAD
from GUD.ORM.tss import TSS

# Initialize session w/ GUD
session = GUDglobals.establish_GUD_session()

# Get all samples in GUD
samples = Sample.select_by_names(
    session,
    treatment=False,
    cell_line=True,
    cancer=True
)

# FANTOM5 enhancer-to-TSS associations
fantom5_file = "/Users/ofornes/Work/mkang/human.associations.hdr.txt"

# TSSs
tsss = TSS.select_all_genic_tss(
    session,
    as_genomic_feature=True
)
# TSS file
tss_file = "tss.bed"
# For each TSS...
for tss in tsss:
    GUDglobals.write(tss_file, tss)

#tss_file = "tss.bed"
## For each line...
#for line in GUDglobals.parse_tsv_file(fantom5_file):
#    for peak in line[1].split(","):
#        m = re.search(
#            "p(\d+)@(\w+)",
#            peak
#        )
#        if m:
#            gene = m.group(2)
#            tss = m.group(1)
#            tss = TSS.select_by_gene_tss(
#                session,
#                m.group(2),
#                m.group(1),
#                as_genomic_feature=True
#            )
#            GUDglobals.write(tss_file, tss)

# For each sample...
for sample in samples:
 
    # Select TADs
    tads = TAD.select_by_sample(
        session,
        sample=sample.name,
        as_genomic_feature=True
    )

    # Select Enhancers
    enhancers = Enhancer.select_by_sample(
        session,
        sample=sample.name,
        as_genomic_feature=True
    )

    if tads and enhancers:

        # Get unique MD5 id
        h = hashlib.md5()
        h.update(sample.name.encode("utf-8"))

        # Create sample directory
        sample_dir = h.hexdigest()
        if not os.path.exists(sample_dir):
            os.makedirs(sample_dir)

        # Sample file
        sample_file = os.path.join(
            sample_dir,
            "sample.txt"
        )
        GUDglobals.write(sample_file, sample.name)

        # TADs file
        tads_file = os.path.join(
            sample_dir,
            "tads.bed"
        )
        # For each TAD...
        for tad in tads:
            GUDglobals.write(tads_file, tad)

        # Enhancers
        enhancers_file = os.path.join(
            sample_dir,
            "enhancers.bed"
        )
        # For each Enhancer...
        for enhancer in enhancers:
            GUDglobals.write(enhancers_file, enhancer)