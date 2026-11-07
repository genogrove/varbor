"""The reference annotation as a pygenogrove graph.

Transcripts are the indexed keys, one index per chromosome, so an interval
query returns the transcripts overlapping it. Exons are external keys: they
take part in the graph but are never returned by a query. Each distinct exon
coordinate is one key, shared by every transcript that uses it.

Edges walk each transcript in genomic order: transcript -> its first exon, and
exon -> next exon. An exon pair is linked once, with the list of transcripts
joining those two exons as the edge's payload, so every edge is one annotated
intron and a transcript is walked by following the edges that list it.

Coordinates are pygenogrove's: 0-based, both ends inclusive.
"""

import collections
import gzip
import re

import pygenogrove as pg

ATTRIBUTE = re.compile(r'(\S+) "([^"]*)"')


def parse_gtf(path):
    """Transcripts of a GTF: {transcript_id: record}, exons sorted genomically.

    A record holds chrom, strand, the 0-based closed span, gene id and name,
    the exons, and for a coding transcript:
      cds_start, cds_frame  the first coding base in transcript orientation
                            (0-based) and that codon's frame
      cds_end               the last coding base in transcript orientation
      start_complete        whether the CDS begins at a full start codon: a
                            start_codon row, or neither a cds_start_NF tag
                            nor a frame other than 0
      stop_codon            whether the annotation gives the stop codon after
                            the CDS (a stop_codon row)

    The file may be gzipped and its lines in any order. A transcript id is
    expected once: a GTF that reuses one at another locus (as copies on the
    pseudoautosomal regions can) would merge both into one record, so it is
    refused, as is an exon whose transcript has no transcript line.
    """
    tx = {}
    exons = collections.defaultdict(list)
    cds = collections.defaultdict(list)
    codons = collections.defaultdict(set)
    opener = gzip.open if str(path).endswith(".gz") else open
    for line in opener(path, "rt"):
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if f[2] not in ("transcript", "exon", "CDS", "start_codon", "stop_codon"):
            continue
        attrs = collections.defaultdict(list)
        for k, v in ATTRIBUTE.findall(f[8]):
            attrs[k].append(v)
        tid = attrs["transcript_id"][0]
        start, end = int(f[3]) - 1, int(f[4]) - 1
        if f[2] == "transcript":
            if tid in tx:
                raise ValueError(f"{path}: transcript {tid} appears more than once")
            tx[tid] = {
                "transcript_id": tid,
                "chrom": f[0],
                "strand": f[6],
                "start": start,
                "end": end,
                "gene_id": attrs["gene_id"][0],
                "gene_name": (attrs.get("gene_name") or [None])[0],
                "cds_start_nf": "cds_start_NF" in attrs.get("tag", []),
            }
        elif f[2] == "exon":
            exons[tid].append((start, end))
        elif f[2] == "CDS":
            cds[tid].append((start, end, int(f[7])))
        else:
            codons[tid].add(f[2])

    orphans = set(exons) - set(tx)
    if orphans:
        raise ValueError(f"{path}: exons of {len(orphans)} transcript(s) without a transcript line, "
                         f"e.g. {sorted(orphans)[0]}")
    for tid, rec in tx.items():
        rec["exons"] = sorted(exons.get(tid, []))
        parts = sorted(cds.get(tid, []))
        if not parts:
            rec["cds_start"] = rec["cds_frame"] = rec["cds_end"] = None
            rec["start_complete"] = rec["stop_codon"] = False
            continue
        if rec["strand"] == "+":
            rec["cds_start"], rec["cds_frame"], rec["cds_end"] = parts[0][0], parts[0][2], parts[-1][1]
        else:
            rec["cds_start"], rec["cds_frame"], rec["cds_end"] = parts[-1][1], parts[-1][2], parts[0][0]
        rec["start_complete"] = ("start_codon" in codons[tid]
                                 or (not rec["cds_start_nf"] and rec["cds_frame"] == 0))
        rec["stop_codon"] = "stop_codon" in codons[tid]
    return tx


def build(gtf_path):
    """Grove of the GTF's transcripts and their shared exons."""
    transcripts = parse_gtf(gtf_path)
    grove = pg.Grove()
    exon_keys = {}
    pairs = collections.defaultdict(list)

    for tid, rec in transcripts.items():
        if not rec["exons"]:
            continue
        payload = {k: v for k, v in rec.items() if k != "exons"}
        tx_key = grove.insert(rec["chrom"],
                              pg.GenomicCoordinate(rec["strand"], rec["start"], rec["end"]),
                              payload)
        chain = []
        for start, end in rec["exons"]:
            ident = (rec["chrom"], rec["strand"], start, end)
            if ident not in exon_keys:
                exon_keys[ident] = grove.add_external_key(
                    pg.GenomicCoordinate(rec["strand"], start, end))
            chain.append(ident)
        grove.add_edge(tx_key, exon_keys[chain[0]], tid)
        for a, b in zip(chain, chain[1:]):
            pairs[(a, b)].append(tid)

    # no call updates an edge's payload, so each pair is linked once, after
    # every transcript using it is known
    for (a, b), tids in pairs.items():
        grove.add_edge(exon_keys[a], exon_keys[b], tids)
    return grove


def exons(grove, tx_key):
    """The transcript's exons in genomic order, as (start, end) pairs."""
    tid = tx_key.data["transcript_id"]
    nxt = grove.get_neighbors(tx_key)
    out = []
    while nxt:
        key = nxt[0]
        out.append((key.value.start, key.value.end))
        nxt = grove.get_neighbors_if(key, lambda tids: tid in tids)
    return out


def transcripts(grove, chrom, strand, start, end):
    """Transcripts overlapping [start, end] on the strand, with their exons.

    Yields (record, exons): the transcript's payload and its exon list.
    """
    for key in grove.intersect(pg.GenomicCoordinate(strand, start, end), chrom):
        yield key.data, exons(grove, key)
