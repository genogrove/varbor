# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- **Nonsense-mediated decay judged on the changed exon structure**: every mutated or truncated effect carries a decay call made in mRNA space on the exons after the event. The call gives the stop's exon, its distance to the exon's 3′ end and the escape rule, so a spliced isoform's moved junctions are taken into account. ([#5](https://github.com/genogrove/varbor/pull/5))

## [0.1.0] - 2026-11-03

### Added

- **Splicing events applied to annotated transcripts**: the annotation is read into a pygenogrove graph of transcripts and shared exons, and SplAdder events of every class are applied with one rule. The isoform an annotated transcript carries is replaced by the one none carries, and the result is translated from the annotated coding start, honouring the CDS phase and the mitochondrial code. ([#1](https://github.com/genogrove/varbor/pull/1))
