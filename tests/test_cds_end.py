"""Mutants of a transcript whose annotation gives no stop codon after the CDS.

Events are built on each transcript's own exons by moving one exon's start:
- ENST00000308647 (+): 3 nt inside the CDS keeps the frame, and 6 intron nt
  kept after a donor bring a stop codon before the CDS end; 1 nt in a 3' UTR
  exon, or 10 nt spliced out of the 5' UTR, changes no protein; an exon ending
  before the CDS end loses it;
- ENST00000530097 (-): moving its second exon's start by 1 nt shifts the frame,
  which then meets its first stop past the annotated CDS end.
"""

import pytest

from varbor import splicing, translate


@pytest.fixture
def no_stop(transcripts, cases):
    r = transcripts[cases["no_stop_codon"]]
    assert r["strand"] == "+" and not r["stop_codon"]
    # an exon inside the CDS, not the first, so its acceptor can move
    k = next(i for i, (s, e) in enumerate(r["exons"])
             if i > 0 and s > r["cds_start"] + 10 and e < r["cds_end"] - 10)
    return r, k


def acceptor_shift(r, k, by):
    source = [r["exons"][k - 1], r["exons"][k]]
    target = [r["exons"][k - 1], (r["exons"][k][0] + by, r["exons"][k][1])]
    return source, target


def test_an_inframe_mutant_ends_at_the_cds_end_with_an_unknown_decay_call(no_stop, fetch):
    r, k = no_stop
    effect = splicing.effect(fetch, r, r["exons"], *acceptor_shift(r, k, 3))
    assert effect.outcome == "mutated" and effect.frameshift is False
    # one codon fewer, and both end at the CDS end
    assert len(effect.mt_protein) == len(effect.wt_protein) - 1
    assert effect.nmd.status == "unknown"


def test_an_inframe_mutant_that_stops_before_the_cds_end_gets_a_decay_call(no_stop, fetch):
    r, _ = no_stop
    k = 4
    before, exon = r["exons"][k - 1], r["exons"][k]
    assert before[0] > r["cds_start"] and exon[1] < r["cds_end"]
    # 6 intron bases kept after the donor bring an in-frame stop codon
    effect = splicing.effect(fetch, r, r["exons"], [before, exon], [(before[0], before[1] + 6), exon])
    assert effect.frameshift is False
    assert len(effect.mt_protein) < len(effect.wt_protein)
    assert effect.nmd.status == "decay"


def test_a_frameshifted_mutant_reads_on_past_the_cds_end(fetch, transcripts, cases):
    r = transcripts[cases["frameshift_read_on"]]
    assert r["strand"] == "-" and not r["stop_codon"]
    effect = splicing.effect(fetch, r, r["exons"], *acceptor_shift(r, 1, 1))
    assert effect.outcome == "mutated" and effect.frameshift is True
    # the new frame runs past the annotated CDS end, to its first stop
    start = effect.cds_offset + r["cds_frame"]
    stop = translate.stop_codon(effect.mrna, start, r["chrom"])
    assert effect.mt_protein == translate.translate(effect.mrna[start:stop], r["chrom"])
    assert effect.nmd.status == "escape"


def test_a_change_in_the_3prime_utr_leaves_the_protein_unchanged(no_stop, fetch):
    r, _ = no_stop
    k = next(i for i, (s, e) in enumerate(r["exons"]) if s > r["cds_end"])
    effect = splicing.effect(fetch, r, r["exons"], *acceptor_shift(r, k, 1))
    assert effect.outcome == "unchanged"


def test_a_changed_start_codon_is_read_literally(fetch, transcripts, cases):
    r = transcripts[cases["split_start_codon"]]
    assert r["strand"] == "-" and r["start_complete"]
    k = next(i for i, (s, e) in enumerate(r["exons"]) if s <= r["cds_start"] <= e)
    exon, before = r["exons"][k], r["exons"][k - 1]
    assert r["cds_start"] - exon[0] < 2
    # the start codon spans this junction; moving the exon's start by 1 nt
    # keeps the CDS start but changes the codon
    effect = splicing.effect(fetch, r, r["exons"], [before, exon], [before, (exon[0] + 1, exon[1])])
    start = effect.cds_offset + r["cds_frame"]
    assert effect.wt_protein[0] == "M"
    assert effect.mt_protein[0] == translate.translate(effect.mrna[start:start + 3], r["chrom"]) != "M"


def test_a_record_from_an_older_graph_is_refused(fetch, transcripts, cases):
    r = dict(transcripts[cases["no_stop_codon"]])
    del r["stop_codon"]
    with pytest.raises(ValueError, match="rebuild"):
        translate.protein(fetch, r, r["exons"])


def test_a_change_in_the_5prime_utr_leaves_the_protein_unchanged(no_stop, fetch):
    r, _ = no_stop
    first = r["exons"][0]
    assert first[0] < r["cds_start"] - 84
    # splicing 10 nt out of the 5' UTR moves the CDS start and end together
    source = [first]
    target = [(first[0], r["cds_start"] - 84), (r["cds_start"] - 73, first[1])]
    effect = splicing.effect(fetch, r, r["exons"], source, target)
    assert effect.outcome == "unchanged"


def test_a_mutant_that_loses_the_cds_end_reads_on_to_its_first_stop(no_stop, fetch):
    r, _ = no_stop
    k = next(i for i, (s, e) in enumerate(r["exons"]) if s <= r["cds_end"] <= e)
    exon, after = r["exons"][k], r["exons"][k + 1]
    # the exon now ends 5 nt before the CDS end
    effect = splicing.effect(fetch, r, r["exons"], [exon, after],
                             [(exon[0], r["cds_end"] - 5), after])
    assert effect.outcome == "mutated"
    start = effect.cds_offset + r["cds_frame"]
    stop = translate.stop_codon(effect.mrna, start, r["chrom"])
    assert stop is not None
    assert effect.mt_protein == translate.translate(effect.mrna[start:stop], r["chrom"])
    assert effect.nmd.status == "decay"
