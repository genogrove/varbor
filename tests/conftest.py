"""Fixtures built from real data.

- annotation.gtf: Ensembl 111 rows (transcript, exon, CDS) of every transcript
  the test events overlap, plus an incomplete-CDS and a mitochondrial transcript.
- spladder/: SplAdder confirmed events (GFF3 and txt) from a TESLA tumor sample.
- genome_slices.fa: GRCh38 bases covering those transcripts, named
  chrom:start-end (0-based, both ends inclusive).
- ensembl_proteins.fa: Ensembl's protein for each fixture transcript.
- cases.json: per event, the transcript and the peptide a junction-based
  reference implementation reports for it.
"""

import json
from pathlib import Path

import pytest

from varbor import graph

DATA = Path(__file__).parent / "data"


def read_fasta(path):
    records, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            records[name] = []
        elif name:
            records[name].append(line)
    return {k: "".join(v) for k, v in records.items()}


@pytest.fixture(scope="session")
def fetch():
    """fetch(chrom, start, end) over the genome slices; refuses what they don't cover."""
    slices = []
    for name, seq in read_fasta(DATA / "genome_slices.fa").items():
        chrom, span = name.rsplit(":", 1)
        start, end = map(int, span.split("-"))
        slices.append((chrom, start, end, seq))

    def get(chrom, start, end):
        for c, s, e, seq in slices:
            if c == chrom and s <= start and end <= e:
                return seq[start - s:end - s + 1]
        raise KeyError(f"{chrom}:{start}-{end} is not in the genome slices")
    return get


@pytest.fixture(scope="session")
def grove():
    return graph.build(DATA / "annotation.gtf")


@pytest.fixture(scope="session")
def transcripts():
    return graph.parse_gtf(DATA / "annotation.gtf")


@pytest.fixture(scope="session")
def ensembl():
    """Ensembl proteins without the leading X of an incomplete start or a trailing stop."""
    return {k: v.lstrip("X").rstrip("*") for k, v in read_fasta(DATA / "ensembl_proteins.fa").items()}


@pytest.fixture(scope="session")
def cases():
    return json.load(open(DATA / "cases.json"))
