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


def genetic_code(chrom):
    return MITOCHONDRIAL if chrom in ("chrMT", "MT", "chrM") else STANDARD


def translate(seq, chrom):
    code = genetic_code(chrom)
    aa = "".join(code.get(seq[i:i + 3], "X") for i in range(0, len(seq) - 2, 3))
    return aa.split("*")[0]


def stop_codon(seq, start, chrom):
    """mRNA index of the first in-frame stop codon from `start`, or None if there is none."""
    code = genetic_code(chrom)
    for i in range(start, len(seq) - 2, 3):
        if code.get(seq[i:i + 3]) == "*":
            return i
    return None


REQUIRED = ("cds_start", "cds_frame", "cds_end", "start_complete", "stop_codon")


def protein(fetch, tx, exons, read_on=False, initiator=None):
    """(protein, mRNA, CDS offset) for the transcript on these exons.

    The protein runs from the annotated CDS start to the first stop codon. It
    is None when the exons no longer hold the CDS start.

    A full start codon is translated as methionine whatever its sequence: the
    initiator tRNA reads CUG or GUG starts as Met. `initiator` is the codon the
    transcript starts with as annotated; where a change has made the first
    codon another one, it is read literally.

    Where the annotation gives no stop codon after the CDS, what the annotated
    frame meets beyond its end is unknown, so translation ends there, unless
    the exons no longer hold that end or `read_on` asks for the first stop
    regardless, as a changed frame does.
    """
    missing = [k for k in REQUIRED if k not in tx]
    if missing:
        raise ValueError(f"transcript {tx.get('transcript_id')} lacks {', '.join(missing)}: "
                         "the annotation graph was built by an older varbor; rebuild it")
    seq = mrna(fetch, tx["chrom"], tx["strand"], exons)
    offset = mrna_offset(tx["strand"], exons, tx["cds_start"])
    if offset is None:
        return None, seq, None
    start = offset + tx["cds_frame"]
    end = cds_end_offset(tx, exons) if not read_on else None
    aa = translate(seq[start:] if end is None else seq[start:end + 1], tx["chrom"])
    if tx["start_complete"] and aa and (initiator is None or seq[start:start + 3] == initiator):
        aa = "M" + aa[1:]
    return aa, seq, offset


def cds_end_offset(tx, exons):
    """mRNA index of the annotated CDS end where translation must stop there, else None.

    That is where the annotation gives no stop codon and the exons still hold
    the CDS end.
    """
    if tx["stop_codon"]:
        return None
    return mrna_offset(tx["strand"], exons, tx["cds_end"])
