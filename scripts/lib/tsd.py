"""Target-site duplication detection, shared by stages 30 and 31.

    empty  : LEFT [tsd] RIGHT
    filled : LEFT [tsd] ELEMENT [tsd] RIGHT

An aligner collapses one copy of the TSD, so the called insert begins or ends
with the same bases that sit beside the empty site.

TSD IS RECORDED, NEVER USED AS A FILTER. Many families make no target-site
duplication at all -- the IS110/IS1111 family is the important case here: it
uses a recombinase-like mechanism, leaves no TSD, and has no terminal inverted
repeats either. Measured on this pipeline's own output: 0/13 IS110-positive
events carry a TSD, against 43/170 for everything else. Gating on TSD would
therefore discard 100% of IS110 events, which is the opposite of the goal.

So a TSD present is positive evidence of transposition; a TSD absent is NOT
evidence against it, and no caller in this pipeline treats it that way.
"""


def find_tsd(ins_seq, left_flank, right_flank, min_tsd=3, max_tsd=30):
    """Longest complexity-filtered direct repeat. Returns (len, seq, side, conf).

    Real TSDs run ~2-9 bp (IS1 and Tn5 make 9, IS3 family 3-4), short enough
    that chance matches are common, so low-complexity repeats are rejected and
    every call carries a confidence.
    """
    if not ins_seq:
        return 0, "", ".", "none"

    def complex_enough(sub):
        if len(set(sub)) < 3:                      # homopolymer / 2-letter
            return False
        return max(sub.count(b) for b in set(sub)) / len(sub) <= 0.75

    hi = min(max_tsd, len(ins_seq))
    for n in range(hi, min_tsd - 1, -1):
        cands = []
        if len(right_flank) >= n and ins_seq[:n].upper() == right_flank[:n].upper():
            cands.append((ins_seq[:n].upper(), "right"))
        if len(left_flank) >= n and ins_seq[-n:].upper() == left_flank[-n:].upper():
            cands.append((ins_seq[-n:].upper(), "left"))
        for sub, side in cands:
            if not complex_enough(sub):
                continue
            return n, sub, side, ("high" if n >= 6 else "medium" if n >= 4 else "low")
    return 0, "", ".", "none"


def diff_alleles(short, long_):
    """Anchor the two allele sequences at both ends; the middle is the insert.

    Returns (insert_seq, left_anchor_len, right_anchor_len).
    """
    if not long_:
        return "", 0, 0
    n = min(len(short), len(long_))
    i = 0
    while i < n and short[i].upper() == long_[i].upper():
        i += 1
    j = 0
    while (j < min(len(short) - i, len(long_) - i)
           and short[len(short) - 1 - j].upper() == long_[len(long_) - 1 - j].upper()):
        j += 1
    return long_[i:len(long_) - j], i, j
