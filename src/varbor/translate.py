"""Transcript sequence and translation.

Sequences come from a `fetch(chrom, start, end)` callable returning the bases
of the 0-based closed interval [start, end]; `genome()` builds one from an
indexed FASTA.
"""

import pygenogrove as pg


def genome(path):
    """A fetch over an indexed FASTA (the index is created if missing)."""
    index = pg.FastaIndex(str(path))
    return lambda chrom, start, end: index.fetch(chrom, start, end + 1)


STANDARD = dict(zip(
    (a + b + c for a in "TCAG" for b in "TCAG" for c in "TCAG"),
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"))
# vertebrate mitochondrial code
MITOCHONDRIAL = dict(STANDARD, TGA="W", AGA="*", AGG="*", ATA="M")
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def mrna(fetch, chrom, strand, exons):
    """The spliced sequence, 5'->3' in transcript orientation."""
    seq = "".join(fetch(chrom, s, e) for s, e in exons).upper()
    return seq[::-1].translate(COMPLEMENT) if strand == "-" else seq


def mrna_offset(strand, exons, pos):
    """The mRNA index of genomic position pos, or None if no exon holds it."""
    ordered = exons if strand == "+" else exons[::-1]
    offset = 0
    for s, e in ordered:
        if s <= pos <= e:
            return offset + (pos - s if strand == "+" else e - pos)
        offset += e - s + 1
    return None


def translate(seq, chrom):
    code = MITOCHONDRIAL if chrom in ("chrMT", "MT", "chrM") else STANDARD
    aa = "".join(code.get(seq[i:i + 3], "X") for i in range(0, len(seq) - 2, 3))
    return aa.split("*")[0]


def stop_codon(seq, start, chrom):
    """mRNA index of the first in-frame stop codon from `start`, or None if there is none."""
    code = MITOCHONDRIAL if chrom in ("chrMT", "MT", "chrM") else STANDARD
    for i in range(start, len(seq) - 2, 3):
        if code.get(seq[i:i + 3]) == "*":
            return i
    return None


def protein(fetch, tx, exons):
    """(protein, mRNA, CDS offset) for the transcript on these exons.

    The protein runs from the annotated CDS start to the first stop codon, read
    on into the 3' UTR if the coding sequence changed. It is None when the
    exons no longer hold the CDS start.
    """
    seq = mrna(fetch, tx["chrom"], tx["strand"], exons)
    offset = mrna_offset(tx["strand"], exons, tx["cds_start"])
    if offset is None:
        return None, seq, None
    return translate(seq[offset + tx["cds_frame"]:], tx["chrom"]), seq, offset


