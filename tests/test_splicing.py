"""Real SplAdder events applied to real transcripts.

Each case is one event per class, direction and strand. Its peptide is the one
a junction-based reference implementation reports for that transcript.
"""

import pytest

from varbor import graph, splicing

from conftest import DATA


@pytest.fixture(scope="module")
def events():
    return {e.event_id: e for e in splicing.read_events(DATA / "spladder", 3)}


def candidates(grove, event):
    lo = min(iso[0][0] for iso in event.isoforms)
    hi = max(iso[-1][1] for iso in event.isoforms)
    return list(graph.transcripts(grove, event.chrom, event.strand, lo, hi))


def test_every_event_has_two_isoforms_and_its_txt_row(events, cases):
    assert {c["event_id"] for c in cases["events"]} <= set(events)
    for event in events.values():
        assert len(event.isoforms) == 2
        assert event.row is not None and event.row["event_id"] == event.event_id


@pytest.mark.parametrize("case_index", range(16))
def test_event_reproduces_the_reference_peptide(case_index, cases, events, grove, fetch, ensembl):
    case = cases["events"][case_index]
    event = events[case["event_id"]]
    cands = candidates(grove, event)
    pair, reason = splicing.source_and_target(event, [ex for _, ex in cands])
    assert reason is None, case["case"]
    source, target = pair

    # the source is carried by an annotated transcript, the target by none
    assert any(splicing.locate(ex, source, target) for _, ex in cands)
    assert not any(splicing.locate(ex, target, source) for _, ex in cands)

    tx, tx_exons = next((t, ex) for t, ex in cands if t["transcript_id"] == case["transcript_id"])
    effect = splicing.effect(fetch, tx, tx_exons, source, target)
    assert effect.outcome == "mutated"
    assert effect.frameshift == case["frameshift"]
    assert effect.wt_protein == ensembl[tx["transcript_id"]]
    # the reference peptide lies at the change: it covers the first changed residue
    at = effect.mt_protein.find(case["reference_peptide"])
    assert at >= 0
    assert at <= effect.var_start < at + len(case["reference_peptide"])


def test_frameshift_is_only_reported_where_the_protein_changes(events, grove, fetch):
    outcomes = set()
    for event in events.values():
        cands = candidates(grove, event)
        pair, _ = splicing.source_and_target(event, [ex for _, ex in cands])
        for tx, tx_exons in cands:
            effect = splicing.effect(fetch, tx, tx_exons, *pair)
            if effect is None:
                continue
            outcomes.add(effect.outcome)
            if effect.outcome in ("mutated", "truncated wt"):
                assert effect.frameshift in (True, False)
            else:
                assert effect.frameshift is None, (event.event_id, tx["transcript_id"], effect.outcome)
    # the fixtures hold both kinds, so neither branch is checked vacuously
    assert "mutated" in outcomes and outcomes - {"mutated", "truncated wt"}


def test_applying_keeps_the_transcripts_outer_ends():
    # an exon skip whose isoforms reach less far than the transcript's exons
    tx = [(100, 200), (300, 400), (500, 600)]
    source = [(150, 200), (300, 400), (500, 550)]
    target = [(150, 200), (500, 550)]
    span = splicing.locate(tx, source, target)
    assert span == (0, 2)
    assert splicing.apply(tx, span, target) == [(100, 200), (500, 600)]


def test_a_transcript_without_the_source_is_not_located():
    tx = [(100, 200), (300, 400), (500, 600)]
    source = [(150, 200), (320, 400), (500, 550)]   # inner exon starts elsewhere
    target = [(150, 200), (500, 550)]
    assert splicing.locate(tx, source, target) is None


def test_a_retained_intron_counts_as_annotated_only_where_an_exon_spans_it():
    spliced = [(100, 200), (300, 400)]
    retained = [(100, 400)]
    spans = [(50, 450)]
    split = [(50, 200), (300, 450)]
    assert splicing.locate(spans, retained, spliced) == (0, 0)
    assert splicing.locate(split, retained, spliced) is None


def event(*isoforms):
    return splicing.Event("e.1", "test", "chr1", "+", "G", tuple(isoforms), None)


def test_an_isoform_counts_as_annotated_whatever_its_outer_exons_reach():
    # alternative 5' site: a uses the donor at 200, b the one at 150
    a = [(100, 200), (300, 400)]
    b = [(100, 150), (300, 400)]
    t1 = [(180, 200), (300, 400), (500, 600)]    # holds a, but starts past b's donor
    t2 = [(50, 150), (300, 400), (500, 600)]     # holds b
    assert splicing.source_and_target(event(a, b), [t1, t2]) == (None, "annotated")
    assert splicing.source_and_target(event(a, b), [t1]) == ((a, b), None)
    # t1 holds a, but its exon can't be rebuilt to b's donor at 150
    assert splicing.locate(t1, a, b) is None


def test_an_event_no_transcript_holds_is_both_novel():
    a = [(100, 200), (300, 400)]
    b = [(100, 250), (300, 400)]
    assert splicing.source_and_target(event(a, b), [[(100, 180), (300, 400)]]) == (None, "both novel")


def test_isoforms_with_different_limits_are_refused():
    a = [(100, 200), (300, 400)]
    b = [(100, 200), (300, 450)]
    assert splicing.source_and_target(event(a, [(100, 150), (300, 450)]), [a]) == (None, "isoform limits differ")


def test_reading_events_from_a_directory_without_them_fails(tmp_path):
    with pytest.raises(FileNotFoundError):
        splicing.read_events(tmp_path, 3)
    with pytest.raises(FileNotFoundError):
        splicing.read_events(DATA / "spladder", 2)
