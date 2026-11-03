"""Alternative-splicing events applied to annotated transcripts.

SplAdder reports every event, whatever its class, as two isoforms with explicit
exon lists. The source is the isoform some annotated transcript carries; the
target is the other one, which no transcript carries. Which is which is decided
against the annotation itself: SplAdder's per-isoform HasNovelJunction flag is
attached to the wrong isoform for about half of the alternative 3'/5' site
events.

A transcript carries the source when the source's inner structure appears in
it consecutively. The inner exon boundaries must match exactly; the outer ends
are only the limits of the event, so the transcript's own outer exon ends are
kept. Applying the event replaces that stretch of the transcript with the
target's exons. This one rule covers exon skipping, alternative 3'/5' sites, intron
retention in either direction, mutually exclusive and multiple exons.

Coordinates are 0-based, both ends inclusive, exons in genomic order.
"""

import collections
import gzip
import os
from pathlib import Path

from varbor.nmd import nmd
from varbor.translate import protein, stop_codon

CLASSES = ("exon_skip", "alt_3prime", "alt_5prime", "intron_retention",
           "mutex_exons", "mult_exon_skip")

Event = collections.namedtuple("Event", "event_id cls chrom strand gene_id isoforms row")
Effect = collections.namedtuple(
    "Effect", "outcome frameshift exons mrna cds_offset wt_protein mt_protein var_start nmd")

def read_events(spladder_dir, confidence):
    """Events of every class from SplAdder's confirmed GFF3 and txt files.

    The GFF3 gives each event's two isoforms; the txt row with the same
    event id is kept whole in `row`, for read support.
    """
    out, found = [], False
    for cls in CLASSES:
        stem = Path(spladder_dir, f"merge_graphs_{cls}_C{confidence}.confirmed")
        gff = Path(f"{stem}.gff3")
        if not gff.exists():
            continue
        found = True
        rows = {}
        txt = Path(f"{stem}.txt.gz")
        if txt.exists():
            with gzip.open(txt, "rt") as fh:
                header = fh.readline().rstrip("\n").split("\t")
                for line in fh:
                    r = dict(zip(header, line.rstrip("\n").split("\t")))
                    rows[r["event_id"]] = r
        out.extend(parse_gff(gff, cls, rows))
    # no file matched: a wrong directory or confidence level must not look
    # like a sample without splicing events
    if not found:
        raise FileNotFoundError(f"{spladder_dir}: no merge_graphs_*_C{confidence}.confirmed.gff3")
    return out


def parse_gff(path, cls, rows):
    genes, isoforms, exons = {}, {}, collections.defaultdict(list)
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
        if f[2] == "gene":
            genes[attrs["ID"]] = (f[0], f[6], attrs.get("GeneName", "").strip('"'))
        elif f[2] == "mRNA":
            isoforms[attrs["ID"]] = attrs["Parent"]
        elif f[2] == "exon":
            exons[attrs["Parent"]].append((int(f[3]) - 1, int(f[4]) - 1))
    by_event = collections.defaultdict(list)
    for iso_id, event_id in sorted(isoforms.items()):
        by_event[event_id].append(sorted(exons[iso_id]))
    for event_id, isos in by_event.items():
        chrom, strand, gene_id = genes[event_id]
        yield Event(event_id, cls, chrom, strand, gene_id, tuple(isos), rows.get(event_id))


def source_and_target(event, tx_exons):
    """(source, target) exon lists, or (None, reason) when the event can't be applied.

    `tx_exons` are the exon lists of the transcripts overlapping the event, on
    its strand. An isoform is annotated when one of them carries it. Counting
    introns alone is not enough: a retained intron has none, yet the annotation
    may lack the exon spanning it.
    The source is the annotated isoform, the target the other one. An event
    whose isoforms are both annotated is a switch between annotated isoforms,
    whose peptides are self.
    """
    if len(event.isoforms) != 2:
        return None, "not two isoforms"
    a, b = event.isoforms
    a_known = any(carries(ex, a, b) for ex in tx_exons)
    b_known = any(carries(ex, b, a) for ex in tx_exons)
    if a_known and b_known:
        return None, "annotated"
    if not (a_known or b_known):
        return None, "both novel"
    source, target = (a, b) if a_known else (b, a)
    if (source[0][0], source[-1][1]) != (target[0][0], target[-1][1]):
        return None, "isoform limits differ"
    return (source, target), None


def carries(tx_exons, isoform, other):
    """(first, last) indices of the transcript exons holding the isoform, or None.

    The isoform's inner structure must appear consecutively: its first exon's
    end, its inner exons and its last exon's start. Its outer ends are only the
    limits of the event. An isoform of one exon has no inner structure; it is
    held by an exon spanning the intron the other isoform splices out.
    """
    if len(isoform) == 1:
        for i, (start, end) in enumerate(tx_exons):
            if start <= other[0][1] and end >= other[-1][0]:
                return i, i
        return None
    for i, (_, end) in enumerate(tx_exons):
        if end != isoform[0][1]:
            continue
        j = i + len(isoform) - 1
        if j < len(tx_exons) and tx_exons[i + 1:j] == isoform[1:-1] and tx_exons[j][0] == isoform[-1][0]:
            return i, j
        return None
    return None


def locate(tx_exons, source, target):
    """(first, last) indices of the transcript exons carrying the source, or None.

    Beyond carrying the source, the transcript's outer exons must reach the
    target's boundaries, so that the exons rebuilt from them are not inverted.
    """
    span = carries(tx_exons, source, target)
    if span is None:
        return None
    i, j = span
    if tx_exons[i][0] > target[0][1] or tx_exons[j][1] < target[-1][0]:
        return None
    return span


def apply(tx_exons, span, target):
    """The transcript's exons with the located stretch replaced by the target's."""
    i, j = span
    inner = list(target)
    inner[0] = (tx_exons[i][0], inner[0][1])
    inner[-1] = (inner[-1][0], tx_exons[j][1])
    return tx_exons[:i] + inner + tx_exons[j + 1:]


def effect(fetch, tx, tx_exons, source, target):
    """The event's effect on one transcript, or None if it doesn't carry the source.

    `frameshift` and `nmd` are set only when the change reaches the protein, as
    a mutated or truncated product; they are None for a noncoding transcript, a
    lost CDS start or a change confined to a UTR. Once it reaches the protein,
    the frame downstream shifts by the whole change in transcript length, and
    decay is judged on the changed exon structure.
    """
    span = locate(tx_exons, source, target)
    if span is None:
        return None
    exons = apply(tx_exons, span, target)
    if tx["cds_start"] is None:
        return Effect("noncoding", None, exons, None, None, None, None, None, None)
    wt, _, _ = protein(fetch, tx, tx_exons)
    mt, seq, offset = protein(fetch, tx, exons)
    if mt is None:
        return Effect("cds start lost", None, exons, seq, None, wt, None, None, None)
    var_start = len(os.path.commonprefix([wt, mt]))
    if mt == wt:
        return Effect("unchanged", None, exons, seq, offset, wt, mt, var_start, None)
    outcome = "truncated wt" if var_start == len(mt) else "mutated"
    frameshift = (sum(e - s + 1 for s, e in exons)
                  - sum(e - s + 1 for s, e in tx_exons)) % 3 != 0
    stop = stop_codon(seq, offset + tx["cds_frame"], tx["chrom"])
    decay = None if stop is None else nmd(tx["strand"], exons, stop)
    return Effect(outcome, frameshift, exons, seq, offset, wt, mt, var_start, decay)
