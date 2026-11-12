#!/usr/bin/env python3
"""Route A (Prodigal meta) and Route B (six-frame) on the same DNA, same table.

THE TWO ROUTES ANSWER DIFFERENT QUESTIONS AND MUST NOT BE MERGED HERE.
  A  predicts genes.  It uses a model of what a gene looks like, so it can be
     wrong about a boundary and silent about a fragment too short or too odd
     for the model.
  B  does not predict anything.  It enumerates every stretch of DNA that COULD
     be read as protein in some frame.  On its own that is not evidence; it
     becomes evidence only when Pfam hits one of the stretches.

SHARED DECISIONS (recorded in routes.json so a rerun cannot drift)
  * TRANSLATION TABLE 11 IN BOTH ROUTES.  Prodigal's meta mode normally picks
    its own bin, and 3 of its 50 bins are table 4.  Letting it do that would
    mean the two routes read the same codons differently, so A is restricted
    to the 47 table-11 bins.  What the UNRESTRICTED meta mode would have
    chosen is still recorded per sequence (`meta_table_unrestricted`): a table-4
    source is a fact about the sequence worth seeing, not something to apply
    silently.
  * NO -c.  `closed=False` keeps genes that run off either end of the
    extraction, which is the normal case for an insertion boundary.
  * PRODIGAL'S OWN MINIMUMS ARE KEPT (min_gene=90 nt, min_edge_gene=60 nt).
    Route B's floor is 20 aa.  The asymmetry is deliberate and is not corrected
    here -- part of what round 1 measures is what A drops for being short.

ROUTE B RULES, FOLLOWING THE BRIEF
  * All 6 frames, split at stop codons, terminal fragments KEPT even though
    they are incomplete.
  * NO start-codon requirement.  Requiring ATG/GTG/TTG would re-lose exactly
    the boundary-truncated CDS this route exists to find.
  * A codon containing anything but ACGT is a HARD SEPARATOR, never an 'X'
    bridging two sides.  Joining across an uncertain codon invents a protein
    that no reading of the DNA supports.
  * Stops are never translated into the protein, and the first codon is NOT
    rewritten to M -- these fragments make no claim to a start.  `first_codon`
    and `has_start_codon` are recorded so a start convention can be applied
    later, explicitly, where there is evidence for one.

COORDINATES.  `nt_start0`/`nt_end0` are 0-based half-open on the FORWARD strand
of the input record, for both routes and both strands.  `frame` is the offset
within the strand's own sequence, so a '-' fragment's frame is read on the
reverse complement.  Any protein here can be regenerated as
    translate( seq[nt_start0:nt_end0] )          when strand == '+'
    translate( revcomp(seq)[L-nt_end0:L-nt_start0] )   when strand == '-'
"""
import argparse, hashlib, json, os, sys, time

import pyrodigal
from Bio.Data import CodonTable

COMP = bytes.maketrans(b"ACGTacgt", b"TGCAtgca")
TABLE = 11


def revcomp(s):
    return s.translate(COMP)[::-1]


def build_table(tid):
    t = CodonTable.unambiguous_dna_by_id[tid]
    aa = dict(t.forward_table)
    for c in t.stop_codons:
        aa[c] = "*"
    return aa, set(t.start_codons)


def fasta_stream(path):
    rid, desc, parts = None, "", []
    with open(path) as f:
        for line in f:
            if line[0] == ">":
                if rid is not None:
                    yield rid, desc, "".join(parts)
                h = line[1:].rstrip("\n")
                rid = h.split()[0] if h.split() else h
                desc = h[len(rid):].strip()
                parts = []
            else:
                parts.append(line.strip())
    if rid is not None:
        yield rid, desc, "".join(parts)


def sixframe(seq, aa_tbl, starts, min_aa):
    """Yield dicts, one per candidate coding fragment, forward-strand coords."""
    L = len(seq)
    rc = revcomp(seq.encode()).decode()
    out = []
    for strand in ("+", "-"):
        s = seq if strand == "+" else rc
        for frame in (0, 1, 2):
            cur = []
            cur_start = None
            pending_left = "edge"
            i = frame
            while i + 3 <= L:
                cod = s[i:i + 3].upper()
                a = aa_tbl.get(cod)
                if a is None or a == "*":
                    bound = "stop" if a == "*" else "ambiguous_codon"
                    if cur and len(cur) >= min_aa:
                        out.append((strand, frame, cur_start, i, "".join(cur),
                                    pending_left, bound, s[cur_start:cur_start + 3].upper()))
                    cur, cur_start = [], None
                    pending_left = bound
                    i += 3
                    continue
                if cur_start is None:
                    cur_start = i
                cur.append(a)
                i += 3
            if cur and len(cur) >= min_aa:
                out.append((strand, frame, cur_start, i, "".join(cur),
                            pending_left, "edge", s[cur_start:cur_start + 3].upper()))
    recs = []
    for strand, frame, a0, a1, prot, lb, rb, fc in out:
        if strand == "+":
            f0, f1 = a0, a1
        else:
            f0, f1 = L - a1, L - a0
        recs.append(dict(strand=strand, frame=frame, nt_start0=f0, nt_end0=f1,
                         protein=prot, left_bound=lb, right_bound=rb,
                         first_codon=fc, has_start_codon=fc in starts))
    return recs


COLS = ["call_id", "route", "seq_key", "nt_start0", "nt_end0", "strand", "frame",
        "length_nt", "length_aa", "translation_table", "partial_begin", "partial_end",
        "start_type", "left_bound", "right_bound", "first_codon", "has_start_codon",
        "score", "confidence", "rbs_motif", "gc_cont", "meta_table_unrestricted"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--out", required=True, help="output prefix (directory must exist)")
    ap.add_argument("--min-aa", type=int, default=20)
    ap.add_argument("--table", type=int, default=TABLE)
    a = ap.parse_args()
    t0 = time.time()

    aa_tbl, starts = build_table(a.table)
    bins = pyrodigal.MetagenomicBins(
        [b for b in pyrodigal.METAGENOMIC_BINS
         if b.training_info.translation_table == a.table])
    gf = pyrodigal.GeneFinder(meta=True, metagenomic_bins=bins, closed=False, mask=False)
    gf_free = pyrodigal.GeneFinder(meta=True, closed=False, mask=False)

    fa = open(a.out + "_prodigal.faa", "w")
    fb = open(a.out + "_sixframe.faa", "w")
    tf = open(a.out + "_calls.tsv", "w")
    tf.write("\t".join(COLS) + "\n")

    n = {"seq": 0, "A": 0, "B": 0, "nontable": 0, "too_short": 0}
    for rid, desc, seq in fasta_stream(a.fasta):
        n["seq"] += 1
        key = rid.split("|", 1)[1] if "|" in rid else rid
        if len(seq) < 60:
            # Prodigal refuses sequences shorter than its edge-gene minimum;
            # route B still runs, so the record is kept, not dropped.
            n["too_short"] += 1
            genes, tbl_free = [], ""
        else:
            genes = gf.find_genes(seq)
            try:
                tbl_free = gf_free.find_genes(seq).training_info.translation_table
            except Exception:
                tbl_free = ""
            if tbl_free and int(tbl_free) != a.table:
                n["nontable"] += 1

        for idx, g in enumerate(genes):
            prot = g.translate(include_stop=False)
            cid = f"A|{key}|{idx:03d}"
            n["A"] += 1
            fa.write(f">{cid} {rid} {g.begin-1}-{g.end} strand={'+' if g.strand>0 else '-'}\n{prot}\n")
            row = dict(call_id=cid, route="prodigal", seq_key=key,
                       nt_start0=g.begin - 1, nt_end0=g.end,
                       strand="+" if g.strand > 0 else "-",
                       frame=(g.begin - 1) % 3 if g.strand > 0 else (len(seq) - g.end) % 3,
                       length_nt=g.end - g.begin + 1, length_aa=len(prot),
                       translation_table=a.table,
                       partial_begin=int(g.partial_begin), partial_end=int(g.partial_end),
                       start_type=g.start_type, left_bound="", right_bound="",
                       first_codon="", has_start_codon="",
                       score=f"{g.score:.3f}", confidence=f"{g.confidence():.2f}",
                       rbs_motif=g.rbs_motif or "", gc_cont=f"{g.gc_cont:.4f}",
                       meta_table_unrestricted=tbl_free)
            tf.write("\t".join(str(row[c]) for c in COLS) + "\n")

        for f in sixframe(seq, aa_tbl, starts, a.min_aa):
            cid = f"B|{key}|{f['strand']}{f['frame']}|{f['nt_start0']}-{f['nt_end0']}"
            n["B"] += 1
            fb.write(f">{cid} {rid}\n{f['protein']}\n")
            row = dict(call_id=cid, route="sixframe", seq_key=key,
                       nt_start0=f["nt_start0"], nt_end0=f["nt_end0"],
                       strand=f["strand"], frame=f["frame"],
                       length_nt=f["nt_end0"] - f["nt_start0"],
                       length_aa=len(f["protein"]), translation_table=a.table,
                       partial_begin=int(f["left_bound"] == "edge"),
                       partial_end=int(f["right_bound"] == "edge"),
                       start_type="", left_bound=f["left_bound"],
                       right_bound=f["right_bound"], first_codon=f["first_codon"],
                       has_start_codon=int(f["has_start_codon"]),
                       score="", confidence="", rbs_motif="", gc_cont="",
                       meta_table_unrestricted=tbl_free)
            tf.write("\t".join(str(row[c]) for c in COLS) + "\n")

        if n["seq"] % 2000 == 0:
            print(f"  {n['seq']} seqs  A={n['A']}  B={n['B']}", flush=True)

    for h in (fa, fb, tf):
        h.close()
    meta = dict(fasta=a.fasta, fasta_sha256=hashlib.sha256(open(a.fasta, "rb").read()).hexdigest(),
                translation_table=a.table, min_aa_sixframe=a.min_aa,
                prodigal_meta_bins=len(bins), prodigal_closed=False,
                prodigal_min_gene=90, prodigal_min_edge_gene=60,
                pyrodigal=pyrodigal.__version__, counts=n,
                seconds=round(time.time() - t0, 1))
    with open(a.out + "_routes.json", "w") as jf:
        json.dump(meta, jf, indent=2)
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
