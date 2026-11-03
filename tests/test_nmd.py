import collections
import re

import pytest

from varbor import graph, splicing, translate
from varbor.nmd import nmd

from conftest import DATA

# three exons of 200, 100 and 300 nt in transcript order, on both strands; on
# the minus strand exon 1 is the rightmost
PLUS = [(1000, 1199), (2000, 2099), (3000, 3299)]


@pytest.mark.parametrize("strand,exons", [("+", PLUS), ("-", [(1000, 1299), (2000, 2099), (3000, 3199)])])
@pytest.mark.parametrize("stop,exon,dist,escape", [
    (99, 1, 100, 2),        # within 100 nt of exon 1's 5' end
    (100, 1, 99, None),     # just past it: exon 1 of 3 is not penultimate
    (249, 2, 50, 3),        # penultimate exon, 50 nt before the last junction
    (248, 2, 51, None),     # 51 nt before it
    (300, 3, 299, 4),       # last exon
])
def test_escape_rules(strand, exons, stop, exon, dist, escape):
    call = nmd(strand, exons, stop)
    status = "decay" if escape is None else "escape"
    assert (call.status, call.exon, call.dist_ejc, call.escape) == (status, exon, dist, escape)


def test_a_single_exon_transcript_escapes():
    assert nmd("+", [(1000, 1999)], 500).escape == 1


def test_a_two_exon_transcript_checks_exon_1_as_penultimate():
    # exon 1 is also the penultimate exon, so a stop near its 3' end escapes by rule 3
    assert nmd("+", [(1000, 1299), (2000, 2299)], 260).escape == 3
    assert nmd("+", [(1000, 1299), (2000, 2299)], 200).escape is None


def test_a_frame_without_a_stop_is_nonstop():
    assert nmd("+", PLUS, None) == ("nonstop", None, None, None)


def test_a_stop_beyond_the_exons_has_no_call():
    assert nmd("+", PLUS, 600) is None


def test_stop_codon_lands_on_ensembls_stop(fetch, transcripts):
    """The base after the last CDS base begins the stop codon (Ensembl's CDS excludes it)."""
    cds_ends, incomplete_end = collections.defaultdict(list), set()
    for line in open(DATA / "annotation.gtf"):
        f = line.split("\t")
        tid = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
        if 'tag "cds_end_NF"' in f[8]:
            incomplete_end.add(tid)
        if f[2] == "CDS":
            cds_ends[tid].append((int(f[3]) - 1, int(f[4]) - 1))
    checked = 0
    for tid, r in transcripts.items():
        if r["cds_start"] is None or tid in incomplete_end:
            continue
        parts = sorted(cds_ends[tid])
        stop_genomic = parts[-1][1] + 1 if r["strand"] == "+" else parts[0][0] - 1
        expected = translate.mrna_offset(r["strand"], r["exons"], stop_genomic)
        if expected is None:          # stop codon split off past the last exon
            continue
        seq = translate.mrna(fetch, r["chrom"], r["strand"], r["exons"])
        offset = translate.mrna_offset(r["strand"], r["exons"], r["cds_start"])
        assert translate.stop_codon(seq, offset + r["cds_frame"], r["chrom"]) == expected, tid
        checked += 1
    assert checked > 50


def genomic_positions(strand, exons):
    """Every exonic genomic position, in transcript order: index i is mRNA position i."""
    positions = [p for s, e in exons for p in range(s, e + 1)]
    return positions if strand == "+" else positions[::-1]


def test_every_case_is_judged_on_its_changed_structure(cases, grove, fetch):
    """The stop is placed in the genome independently, base by base, and the
    changed exon holding it must be the one the call names."""
    events = {e.event_id: e for e in splicing.read_events(DATA / "spladder", 3)}
    judged = 0
    for case in cases["events"]:
        event = events[case["event_id"]]
        lo = min(iso[0][0] for iso in event.isoforms)
        hi = max(iso[-1][1] for iso in event.isoforms)
        cands = list(graph.transcripts(grove, event.chrom, event.strand, lo, hi))
        (source, target), _ = splicing.source_and_target(event, [ex for _, ex in cands])
        tx, ex = next((t, e) for t, e in cands if t["transcript_id"] == case["transcript_id"])
        effect = splicing.effect(fetch, tx, ex, source, target)
        stop = translate.stop_codon(effect.mrna, effect.cds_offset + tx["cds_frame"], tx["chrom"])
        if stop is None:
            assert effect.nmd.status == "nonstop", case["case"]
            continue
        assert effect.nmd.status in ("decay", "escape"), case["case"]
        where = genomic_positions(tx["strand"], effect.exons)[stop]
        ordered = effect.exons if tx["strand"] == "+" else effect.exons[::-1]
        number, (s, e) = next((n, x) for n, x in enumerate(ordered, 1) if x[0] <= where <= x[1])
        assert effect.nmd.exon == number, case["case"]
        assert effect.nmd.dist_ejc == (e - where if tx["strand"] == "+" else where - s), case["case"]
        judged += 1
    assert judged >= 12
