# varbor

Variant effect prediction using [genogrove](https://github.com/genogrove).

varbor builds the reference annotation as a pygenogrove graph. Transcripts are the indexed keys, and exons are shared between transcripts and linked in transcript order. A variant is applied to the transcripts that carry it, as a change to the transcript's structure or its sequence. The result is translated from the annotated coding start.

Planned variant sources, in order:
1. alternative splicing (SplAdder events);
2. exitrons;
3. gene fusions;
4. SNVs and indels.

## Status

Early development: releases are `0.x`, and the API may change between them. Alternative splicing is supported; the other variant sources follow.

## License

GPL-3.0, see [LICENSE](LICENSE).
