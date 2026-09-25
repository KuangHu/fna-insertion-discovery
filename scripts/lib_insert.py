"""One place to get the inserted sequence out of a locus row.

WHY THIS EXISTS -- bug found 2026-09-24 by code review of 0262ce3.

`lcp_bp` is an offset in the SHORT allele.  That is its correct and intended
meaning: it is the insertion point in the empty target, the offset that makes

    target_seq[:offset] + insert_seq + target_seq[offset:] == derived allele

true.  86_catalogue.py, 99_mobility_join.py and 100_event_dedup.py each
independently did `long_allele[lcp : lcp + inserted_len]`, using a short-allele
offset as a long-allele coordinate.

On the EXACT path the two coincide -- lcp is a common prefix, so it counts the
same bases in both alleles -- and the re-slice is right.  On the TOLERANT path
tolerant_decompose() returns `"lcp": sp` (short frame) while the sequence it
hashed into `inserted_md5` was `long_[bl0:bl1]` (long frame).  sp and bl0
differ by exactly (upstream insertions - upstream deletions) in the alignment,
so any indel before the junction shifts the re-slice by that many bases.

Measured on 19 species, 199,766 catalogue rows:
    exact               125,947 rows        0 wrong   0.00%
    tolerant_alignment   73,819 rows   22,972 wrong  31.12%
Synthetic dose-response confirms the mechanism: a d bp upstream deletion gives
shift -d, a d bp upstream insertion gives +d, and the quality gate passes both.

`inserted_md5` was hashed from the CORRECT sequence, so it is ground truth even
where the re-slice is wrong.  That makes the damage repairable from stored data
and Layer 2 does not have to be re-run.  Post-fix runs of 70_ also write
`insert_start_in_long_bp` explicitly and skip the search.
"""
import collections
import hashlib
import os


def md5(s):
    return hashlib.md5(s.encode()).hexdigest()[:12] if s else "EMPTY"


def rc(s):
    return s.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def canon_key(s):
    """Orientation-invariant key for a sequence.

    `inserted_md5` is a plain hash of the allele as stored, so the same element
    inserted on opposite strands hashes to two different values.  Anything that
    asks "is this the same element somewhere else" must use this instead.
    """
    if not s:
        return ""
    r = rc(s)
    return hashlib.md5((s if s <= r else r).encode()).hexdigest()[:16]


def canonical_insert_key(long_seq, start, ilen):
    """Orientation- AND junction-invariant key for an insert.

    WHY canon_key IS NOT ENOUGH -- review finding #2, confirmed 2026-09-24.

    When the junction carries a direct repeat of length d, the same physical
    insertion can be written d+1 different ways: the boundary may sit anywhere
    inside the repeat and every choice reconstructs the identical derived
    allele. The alignment picks one, and which one it picks is not symmetric
    under reverse-complementation. Measured on 2,961 E. coli loci, re-choosing
    the junction in the revcomp frame:

        junction overlap 0      79.7% of loci keep the same key
        junction overlap 1-3     0.0%
        junction overlap 4-15    0.0%
        junction overlap >15     0.0%

    So any locus with microhomology -- which is most of them, and all the ones
    carrying a TSD -- keys differently depending on which strand its assembly
    happened to be deposited on. Consequence in the shipped output: 351 of
    21,051 E. coli events (1.67%) are the same event counted twice, by dedup's
    own definition of an event.

    The fix is lossless rather than lossy: slide the insert across the whole
    direct repeat, and take the smallest canonical key over the equivalent
    placements. Trimming the ends would also work but would throw away the
    ability to tell apart elements that differ only at their termini.
    """
    if ilen <= 0 or not long_seq:
        return ""
    lo = hi = start
    while lo > 0 and long_seq[lo - 1] == long_seq[lo - 1 + ilen]:
        lo -= 1
    while hi + ilen < len(long_seq) and long_seq[hi] == long_seq[hi + ilen]:
        hi += 1
    return min(canon_key(long_seq[s:s + ilen]) for s in range(lo, hi + 1))


def read_alleles(path):
    """Parse an alleles.fna written by 70_: >LOCUS|ALLELE|len=..|n=.."""
    out, cur = collections.defaultdict(dict), None
    if not os.path.exists(path):
        return out
    for line in open(path):
        if line.startswith(">"):
            p = line[1:].strip().split("|")
            cur = (p[0], p[1]) if len(p) >= 2 else None
            if cur:
                out[cur[0]][cur[1]] = []
        elif cur:
            out[cur[0]][cur[1]].append(line.strip())
    return {k: {a: "".join(v).replace("-", "").upper() for a, v in d.items()}
            for k, d in out.items()}


def extract_insert(row, long_seq, window=None):
    """(seq, status). See locate_insert for the status values."""
    st, status = locate_insert(row, long_seq, window)
    if st < 0:
        return "", status
    return long_seq[st:st + int(row["inserted_len"])], status


def locate_insert(row, long_seq, window=None):
    """Return (long_allele_start, status) for one locus row.

    status is one of
      coord        taken from insert_start_in_long_bp  (post-fix 70_ output)
      lcp          lcp_bp was already the right coordinate (the exact path)
      repaired     lcp_bp was wrong; the offset matching inserted_md5 was found
      unresolved   no offset in the allele reproduces inserted_md5

    The search is nearest-offset-first and by default covers every valid start
    in the long allele. A +/-300 bp window left 164 of 3,676 E. coli repairs
    unresolved -- those are loci with a small lcp and a large insert, where the
    true shift exceeds 300 -- and the scan is bounded by the allele length
    anyway. Measured on 41,607 E. coli loci: 0 cases where more than one offset
    matched, so the md5 anchor is decisive and "nearest" never has to break a
    real tie.

    An `unresolved` row must NOT be exported as sequence.  Emitting a wrong
    insert is worse than emitting none: it silently contaminates CDS calling,
    clustering and every downstream bag.
    """
    try:
        ilen = int(row.get("inserted_len") or 0)
    except ValueError:
        return -1, "unresolved"
    if ilen <= 0 or not long_seq:
        return -1, "unresolved"
    stored = (row.get("inserted_md5") or "").strip()

    explicit = (row.get("insert_start_in_long_bp") or "").strip()
    if explicit not in ("", ".", "NA"):
        try:
            st = int(explicit)
        except ValueError:
            st = None
        if st is not None and 0 <= st and st + ilen <= len(long_seq):
            if not stored or md5(long_seq[st:st + ilen]) == stored:
                return st, "coord"

    try:
        lcp = int(row.get("lcp_bp") or 0)
    except ValueError:
        return -1, "unresolved"
    if 0 <= lcp and lcp + ilen <= len(long_seq):
        if not stored or md5(long_seq[lcp:lcp + ilen]) == stored:
            return lcp, "lcp"

    # md5-anchored repair. Nearest offset first, so the result is deterministic
    # and is the minimal shift rather than an arbitrary match in a repeat.
    if stored:
        hi = len(long_seq) - ilen           # last valid start
        span = window if window is not None else max(lcp, hi - lcp)
        for d in range(1, span + 1):
            for st in (lcp - d, lcp + d):
                if 0 <= st <= hi and md5(long_seq[st:st + ilen]) == stored:
                    return st, "repaired"
    return -1, "unresolved"
