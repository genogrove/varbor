"""Nonsense-mediated decay, judged from where a transcript's stop codon falls.

Positions are in mRNA space, on whatever exon structure the transcript has after
a variant is applied, so a changed structure moves the exon junctions the
decision depends on. A stop escapes decay when no exon junction lies far enough
downstream of it to leave an exon junction complex behind:

  1  the transcript has a single exon;
  2  the stop lies within 100 nt of the 5' end of exon 1, where translation can
     re-initiate downstream;
  3  the stop lies in the penultimate exon, within 50 nt of the last junction;
  4  the stop lies in the last exon.

Otherwise the transcript is a decay target.
"""

import collections

NMD = collections.namedtuple("NMD", "exon dist_ejc escape triggers")


def nmd(strand, exons, stop):
    """Decay call for a stop codon starting at mRNA index `stop`, or None if no exon holds it.

    `exon` is the 1-based exon number in transcript order, `dist_ejc` the number of
    bases from the stop to its exon's 3' end, `escape` the rule that lets it
    escape (None if none does), and `triggers` whether decay is expected.
    """
    lengths = [e - s + 1 for s, e in (exons if strand == "+" else exons[::-1])]
    start = 0
    for number, length in enumerate(lengths, 1):
        if start <= stop < start + length:
            dist_ejc = start + length - 1 - stop
            break
        start += length
    else:
        return None

    last = len(lengths)
    escape = None
    if number == 1:
        if last == 1:
            escape = 1
        elif stop < 100:
            escape = 2
    if escape is None:
        if number == last - 1:
            if dist_ejc <= 50:
                escape = 3
        elif number == last:
            escape = 4
    return NMD(number, dist_ejc, escape, escape is None)
