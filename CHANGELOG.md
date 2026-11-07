# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.2.1] - 2026-11-06

### Fixed

- **Start and end translated the way the annotation defines them**: a full start codon reads as `M` whatever its sequence, unless an event changes it, and without a `stop_codon` row translation ends at the annotated CDS end. A mutant that reaches that end in frame without a stop gets the new decay status `unknown`, because its stop is not known; a mutant whose frame shifts reads on to its first stop. Graphs built by v0.2.0 lack the new fields and must be rebuilt. ([#7](https://github.com/genogrove/varbor/pull/7))

## [0.2.0] - 2026-11-03

### Added

- **Nonsense-mediated decay judged on the changed exon structure**: every mutated or truncated effect carries a decay call made in mRNA space on the exons after the event, so a spliced isoform's moved junctions are taken into account. Its status is `decay`, `escape` (with the rule) or `nonstop`, the last for a reading frame that meets no stop codon. ([#5](https://github.com/genogrove/varbor/pull/5))

## [0.1.0] - 2026-11-03

### Added

- **Splicing events applied to annotated transcripts**: the annotation is read into a pygenogrove graph of transcripts and shared exons, and SplAdder events of every class are applied with one rule. The isoform an annotated transcript carries is replaced by the one none carries, and the result is translated from the annotated coding start, honouring the CDS phase and the mitochondrial code. ([#1](https://github.com/genogrove/varbor/pull/1))
