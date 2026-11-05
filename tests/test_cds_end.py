"""Mutants of a transcript whose annotation gives no stop codon after the CDS.

Events are built on each transcript's own exons by moving one exon's start:
- ENST00000308647 (+): 3 nt inside the CDS keeps the frame; 1 nt in a 3' UTR
  exon changes no protein;
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
    r = transcripts[cases["non_aug_start"]]
    as_annotated, seq, offset = translate.protein(fetch, r, r["exons"])
    codon = seq[offset + r["cds_frame"]:offset + r["cds_frame"] + 3]
    assert as_annotated[0] == "M"
    # the same exons, told the annotated start codon was another one
    changed, _, _ = translate.protein(fetch, r, r["exons"], initiator="ATG")
    assert changed[0] == translate.translate(codon, r["chrom"]) != "M"


def test_a_record_from_an_older_graph_is_refused(fetch, transcripts, cases):
    r = dict(transcripts[cases["no_stop_codon"]])
    del r["stop_codon"]
    with pytest.raises(ValueError, match="rebuild"):
        translate.protein(fetch, r, r["exons"])
