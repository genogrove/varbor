from varbor import translate


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
