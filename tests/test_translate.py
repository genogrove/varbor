from pathlib import Path

from varbor import graph, translate


def test_reference_proteins_match_ensembl(fetch, transcripts, ensembl):
    coding = [t for t, r in transcripts.items() if r["cds_start"] is not None and t in ensembl]
    assert coding
    for tid in coding:
        r = transcripts[tid]
        protein, _, _ = translate.protein(fetch, r, r["exons"])
        assert protein == ensembl[tid], tid


def test_incomplete_cds_start_reads_in_the_phase_the_gtf_gives(fetch, transcripts, ensembl, cases):
    tid = cases["cds_start_nf"]
    r = transcripts[tid]
    assert r["cds_start_nf"] and r["cds_frame"] != 0
    protein, seq, offset = translate.protein(fetch, r, r["exons"])
    assert protein == ensembl[tid]
    # the first base of the CDS is not the first base of a codon here
    assert translate.translate(seq[offset:], r["chrom"]) != ensembl[tid]


def test_mitochondrial_transcript_uses_the_mitochondrial_code(fetch, transcripts, ensembl, cases):
    tid = cases["mitochondrial"]
    r = transcripts[tid]
    assert r["chrom"] == "chrMT"
    protein, seq, offset = translate.protein(fetch, r, r["exons"])
    assert protein == ensembl[tid]


def test_genome_fetch_is_closed_and_zero_based(tmp_path):
    fasta = tmp_path / "g.fa"
    fasta.write_text(">chr1\nACGTACGTAC\n")
    fetch = translate.genome(fasta)
    assert fetch("chr1", 0, 0) == "A"
    assert fetch("chr1", 2, 5) == "GTAC"


def test_a_full_non_aug_start_reads_as_methionine(fetch, transcripts, ensembl, cases):
    tid = cases["non_aug_start"]
    r = transcripts[tid]
    seq = translate.mrna(fetch, r["chrom"], r["strand"], r["exons"])
    offset = translate.mrna_offset(r["strand"], r["exons"], r["cds_start"])
    codon = seq[offset:offset + 3]
    assert r["start_complete"] and codon != "ATG"
    protein, _, _ = translate.protein(fetch, r, r["exons"])
    assert protein == ensembl[tid]
    # read literally, the codon gives another residue
    assert translate.translate(codon, r["chrom"]) != "M"


def test_a_non_aug_start_without_a_start_codon_row_reads_as_methionine(fetch, transcripts, ensembl, cases):
    tid = cases["no_start_codon_row"]
    r = transcripts[tid]
    gtf = (Path(__file__).parent / "data" / "annotation.gtf").read_text().splitlines()
    assert not any("\tstart_codon\t" in l and f'"{tid}"' in l for l in gtf)
    # nothing marks the start as incomplete: no cds_start_NF tag, phase 0
    assert not r["cds_start_nf"] and r["cds_frame"] == 0
    protein, seq, offset = translate.protein(fetch, r, r["exons"])
    assert seq[offset:offset + 3] != "ATG"
    assert protein == ensembl[tid] and protein[0] == "M"


def test_a_nonzero_phase_alone_marks_the_start_incomplete(transcripts, cases, tmp_path):
    tid = cases["cds_start_nf"]
    assert transcripts[tid]["cds_frame"] != 0
    # Ensembl tags every such start; without the tag the phase still tells
    gtf = (Path(__file__).parent / "data" / "annotation.gtf").read_text().splitlines(keepends=True)
    rows = [l.replace(' tag "cds_start_NF";', "") for l in gtf if f'"{tid}"' in l]
    assert not any("cds_start_NF" in l or "\tstart_codon\t" in l for l in rows)
    (tmp_path / "a.gtf").write_text("".join(rows))
    assert graph.parse_gtf(tmp_path / "a.gtf")[tid]["start_complete"] is False


def test_a_cds_without_a_stop_codon_ends_where_the_annotation_does(fetch, transcripts, ensembl, cases):
    tid = cases["no_stop_codon"]
    r = transcripts[tid]
    assert not r["stop_codon"]
    protein, seq, offset = translate.protein(fetch, r, r["exons"])
    assert protein == ensembl[tid]
    # read on to the first stop, it runs into sequence the annotation doesn't vouch for
    assert len(translate.translate(seq[offset + r["cds_frame"]:], r["chrom"])) > len(protein)
