import numpy as np
from Bio.SeqFeature import FeatureLocation, SeqFeature

class GenomicFeature(SeqFeature):
    """
    Implements a Genomic Feature object based on the Biopython's Sequence
    Feature object.

    Attributes:
    chrom {str} chromosome of the feature location {FeatureLocation} location
    of the feature on the genome type {str} the specified type of feature (e.g.
    gene, TSS, repeat...) strand {int} on the DNA sequence. \"1\" indicates the
    plus strand; \"-1\" the minus strand; \"0\" for unknown or not applicable
    strand id {str} identifier for the feature qualifiers {dict} qualifiers of
    the feature profile {array} of scores per nucleotide.
    """

    def __init__(
        self,
        chrom,
        start,
        end,
        score=0,
        strand=None,
        feat_type="Feature",
        feat_id="NA",
        qualifiers=None,
        profile=None
    ):

        self.chrom = chrom
        self._start = start
        self._end = end
        self.location = FeatureLocation(
            self.start,
            self.end
        )
        self.score = score
        self._strand = strand
        self.strand = self.strand_binary
        self.type = feat_type
        self.id = feat_id
        self.qualifiers = qualifiers

        if profile is not None:
            if not isinstance(profile, np.ndarray):
                raise ValueError(
                    "Input profile is not an NumPy array!"
                )

        self.profile = profile

    @property
    def start(self):

        return int(self._start)

    @property
    def start_1_based(self):

        return self.start + 1

    @property
    def end(self):

        return int(self._end)

    @property
    def end_1_based(self):

        return self.end

    @property
    def strand_binary(self):

        if self._strand == "+":
            return 1
        if self._strand == "-":
            return -1

        return 0

    @property
    def strand_string(self):

        if self.strand == 1:
            return "+"
        if self.strand == -1:
            return "-"

        return "."

    @property
    def gud_strand(self):

        return self._strand

    def overlaps(self, feat):
        """
        Returns if feat overlaps this one or not.
        """

        try:
            # Note that, for 0-based feature, if
            # feat1.end == feat2.start, feats don't
            # overlap but are book-ended. E.g.
            #
            # ---------1111111---------
            # -------22222------------- => True
            # ----------22222---------- => True
            # -------------22222------- => True
            # -----22222--------------- => False
            # ---------------22222----- => False
            if self.chrom == feat.chrom and \
                    self.start < feat.end and \
                    self.end > feat.start:

                return True
        except:
            raise ValueError(
                "Could not calculate overlap!"
            )

    def __str__(self):

        return "{}\t{}\t{}\t{}\t{}\t{}".\
            format(
                self.chrom,
                self.start,  # 0-based for BED format
                self.end,
                self.id,
                self.score,
                self.strand_string
            )

    def __repr__(self):

        return "<%s(%s, %s, %s, %s, %s, %s)>" % \
            (
                self.type,
                "chrom={}".format(self.chrom),
                "start={}".format(self.start),
                "end={}".format(self.end),
                "id={}".format(self.id),
                "score={}".format(self.score),
                "strand={}".format(self.strand)
            )

    def serialize(self):
        return {
            "chrom": self.chrom,
            "start": self.start,
            "end": self.end,
            "type": self.type,
            "id": self.id,
            "score": self.score,
            "strand": self.strand,
            "qualifiers": self.qualifiers,
        }
