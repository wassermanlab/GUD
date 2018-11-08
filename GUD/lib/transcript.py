#
# Add compute UTRs (exons - coding exons for 1st and last exons)
# Compute introns
# Compute promoter
#

def compute_exons(self):
    """
    """

    if not self._exons:
        # Initialize
        exons = []
        strand = None
        exon_starts = self.exonStarts.split(",")
        exon_ends = self.exonEnds.split(",")
        if self.strand is None or self.strand == "+" or self.strand == 1:
            strand = 1
        elif self.strand == "-" or self.strand == -1:
            strand = -1
        # For each exon...
        for s, e in zip(exon_starts, exon_ends):
            # Last element of exon_starts and exon_ends
            # is an empty string
            if s != "" and e != "":
                # Initialize exon
                exon_num = len(exons) + 1
                # Add exon to exons
                exons.append(
                    Region(
                        self.chrom,
                        FeatureLocation(int(s), int(e)),
                        id = "Exon{}-{}".format(exon_num, self.name2),
                        type = "exon",
                        strand=strand,
                        qualifiers = {
                            "source" : "ref_gene",
                        }
                    )
                )
        # Add exons
        self._exons = exons

    return self._exons

def compute_coding_exons(self):
    """
    """

    if not self._coding_exons:
        # Initialize
        exons = []
        strand = None
        exon_starts = self.exonStarts.split(",")
        exon_ends = self.exonEnds.split(",")
        if self.strand is None or self.strand == "+" or self.strand == 1:
            strand = 1
        elif self.strand == "-" or self.strand == -1:
            strand = -1
        # For each exon...
        for s, e in zip(exon_starts, exon_ends):
            # Last element of exon_starts and exon_ends
            # is an empty string
            if s != "" and e != "":
                s = int(s)
                e = int(e)
                # Initialize exon
                exon_num = len(exons) + 1
                # If coding exon...
                if s < self.cdsStart:
                    if e > self.cdsStart:
                        if e <= self.cdsEnd:
                            # Add exon to exons
                            exons.append(
                                Region(
                                    self.chrom,
                                    FeatureLocation(self.cdsStart, e),
                                    id = "Exon{}-{}".format(exon_num, self.name2),
                                    type = "exon",
                                    strand=strand,
                                    qualifiers = {
                                        "source" : "ref_gene",
                                    }
                                )
                            )
                        else:
                            exons.append(
                                Region(
                                    self.chrom,
                                    FeatureLocation(self.cdsStart, self.cdsEnd),
                                    id = "Exon{}-{}".format(exon_num, self.name2),
                                    type = "exon",
                                    strand=strand,
                                    qualifiers = {
                                        "source" : "ref_gene",
                                    }
                                )
                            )
                elif s < self.cdsEnd:
                    if e <= self.cdsEnd:
                        exons.append(
                            Region(
                                self.chrom,
                                FeatureLocation(s, e),
                                id = "Exon{}-{}".format(exon_num, self.name2),
                                type = "exon",
                                strand=strand,
                                qualifiers = {
                                    "source" : "ref_gene",
                                }
                            )
                        )
                    else:
                        exons.append(
                            Region(
                                self.chrom,
                                FeatureLocation(s, self.cdsEnd),
                                id = "Exon{}-{}".format(exon_num, self.name2),
                                type = "exon",
                                strand=strand,
                                qualifiers = {
                                    "source" : "ref_gene",
                                }
                            )
                        )
        # Add exons
        self._coding_exons = exons

    return self._coding_exons