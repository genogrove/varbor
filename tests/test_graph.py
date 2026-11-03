import gzip
import random

import pygenogrove as pg
import pytest

from varbor import graph

from conftest import DATA


def test_every_transcript_walks_its_gtf_exons(grove, transcripts):
    for key in grove.intersect(pg.GenomicCoordinate("*", 0, 300_000_000)):
        tid = key.data["transcript_id"]
        assert graph.exons(grove, key) == transcripts[tid]["exons"], tid


def test_exons_are_shared_between_transcripts(grove, transcripts):
    distinct = {(r["chrom"], r["strand"], s, e) for r in transcripts.values() for s, e in r["exons"]}
    rows = sum(len(r["exons"]) for r in transcripts.values())
    assert grove.external_vertex_count() == len(distinct)
    assert len(distinct) < rows


def test_an_exon_pair_is_one_edge_listing_its_transcripts(grove, transcripts):
    pairs = {}
    for tid, r in transcripts.items():
        for a, b in zip(r["exons"], r["exons"][1:]):
            pairs.setdefault((r["chrom"], r["strand"], a, b), set()).add(tid)
    first_exon_edges = len(transcripts)
    assert grove.edge_count() == len(pairs) + first_exon_edges


def test_query_is_strand_aware(grove, transcripts):
    r = next(r for r in transcripts.values() if r["strand"] == "+")
    plus = {t["transcript_id"] for t, _ in graph.transcripts(grove, r["chrom"], "+", r["start"], r["end"])}
    minus = {t["transcript_id"] for t, _ in graph.transcripts(grove, r["chrom"], "-", r["start"], r["end"])}
    assert r["transcript_id"] in plus
    assert r["transcript_id"] not in minus


def test_serialized_graph_answers_the_same(grove, transcripts, tmp_path):
    path = tmp_path / "annotation.gg"
    grove.serialize(str(path))
    view = pg.GroveView.open(str(path))
    for tid, r in transcripts.items():
        hits = [k for k in view.intersect(pg.GenomicCoordinate(r["strand"], r["start"], r["start"]), r["chrom"])
                if k.data["transcript_id"] == tid]
        assert len(hits) == 1, tid
        walked, nxt = [], view.get_neighbors(hits[0])
        while nxt:
            walked.append((nxt[0].value.start, nxt[0].value.end))
            nxt = view.get_neighbors_if(nxt[0], lambda tids: tid in tids)
        assert walked == r["exons"], tid


def test_gtf_may_be_gzipped_shuffled_and_have_blank_lines(transcripts, tmp_path):
    lines = open(DATA / "annotation.gtf").read().splitlines(keepends=True)
    random.Random(1).shuffle(lines)
    path = tmp_path / "shuffled.gtf.gz"
    with gzip.open(path, "wt") as out:
        out.write("\n".join(lines[:10]) + "\n\n" + "".join(lines[10:]))
    assert graph.parse_gtf(path) == transcripts


def test_a_transcript_id_used_twice_is_refused(tmp_path):
    lines = open(DATA / "annotation.gtf").read().splitlines(keepends=True)
    first = next(l for l in lines if l.split("\t")[2] == "transcript")
    path = tmp_path / "twice.gtf"
    path.write_text("".join(lines) + first.replace(first.split("\t")[0], "chrY", 1))
    with pytest.raises(ValueError, match="more than once"):
        graph.parse_gtf(path)


def test_an_exon_without_its_transcript_line_is_refused(tmp_path):
    lines = open(DATA / "annotation.gtf").read().splitlines(keepends=True)
    tid = next(l for l in lines if l.split("\t")[2] == "transcript").split('transcript_id "')[1].split('"')[0]
    path = tmp_path / "orphan.gtf"
    path.write_text("".join(l for l in lines if not (l.split("\t")[2] == "transcript" and f'"{tid}"' in l)))
    with pytest.raises(ValueError, match="without a transcript line"):
        graph.parse_gtf(path)
