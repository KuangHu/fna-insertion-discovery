#!/usr/bin/env python3
"""Decompose the NON-EXTREME allele pairs 70_ never looks at.

70_ decomposes the longest allele against the shortest. Where a locus holds an
intermediate allele, that single comparison describes a compound difference:

    A  L-----------R       A->C reports  X+Y  as one insert
    B  L---X-------R       A->B reports  X
    C  L---X---Y---R       B->C reports  Y

49.9% of loci (127,602 of 255,770) carry at least one intermediate LENGTH
class, so the structure is common. Whether those pairs decompose is a separate
question, and this answers it.

THE CHAIN. Nothing is counted until it survives the previous stage, because
each stage removes a different way of being wrong:

  1 candidate pairs    ordered allele pairs (shorter, longer) at one anchored
                       locus, EXCLUDING the (shortest, longest) pair 70_
                       already did. Pairs of equal length are skipped: they
                       differ by substitution and cannot stack.

  2 decomposed         the pair passes the SAME gates as the extremes --
                       --min-insert 500, coverage, identity, dominance. Not
                       relaxed: a pair admitted under looser gates is not
                       comparable to what is already in the catalogue and would
                       make stage 3 uninterpretable.

  3 novel pairs        after canonical keying -- (canonical insert, canonical
                       +/-100bp context) -- the pair is NOT already an event in
                       the existing catalogue, and not a duplicate of another
                       new pair. This is where a compound difference expressed
                       two ways (A->C vs A->B + B->C) collapses to one event,
                       and it is the only number that means "deeper mining".

TWO CLASSES OF PAIR, kept in separate files because they are different objects:

  empty_target   the shorter allele of the pair IS the locus's shortest allele,
                 so the target really is the pre-insertion site. Directly
                 comparable to a catalogue event, and the only class that may
                 enter the database or a bag flank.

  nested_target  the shorter allele is an INTERMEDIATE allele, so the target
                 ALREADY CARRIES an element. Structurally real -- this is the
                 B->C comparison that recovers Y -- but it is not an
                 (empty site, insert) pair. Measured on the S. pneumoniae
                 pilot: 2,088 of 3,315 novel pairs (63.0%) are this class.
                 Feeding them into the bags would put element sequence into a
                 flank that is defined as pre-insertion.

What this does NOT claim: a nested pair is a STRUCTURAL relation between
alleles at one site. It is not evidence of two independent transposition
events, and the parent/child ordering here is by length, not by ancestry.

PARALLELISM: stage 1-2 is a Biopython global alignment per pair and dominates
the runtime -- 5h21m single-threaded on the SMALLEST species. Panels are
independent, so they run in a process pool; stage 3 stays in the parent because
its dedup is global.
"""
import argparse, collections, csv, glob, hashlib, importlib.util, os, sys
csv.field_size_limit(sys.maxsize)
D = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, D)
import lib_insert as LI
_s = importlib.util.spec_from_file_location(
    "l2", os.path.join(D, "70_allele_reconstructor.py"))
L2 = importlib.util.module_from_spec(_s); _s.loader.exec_module(L2)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--recon", nargs="+", required=True)
    ap.add_argument("--events", required=True, help="existing dedup loci_to_event")
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-insert", type=int, default=500)
    ap.add_argument("--context", type=int, default=100)
    ap.add_argument("--max-pairs-per-locus", type=int, default=12)
    ap.add_argument("--threads", type=int, default=1)
    a = ap.parse_args()

    dirs = []
    for pat in a.recon:
        dirs.extend(sorted(glob.glob(pat)) if any(c in pat for c in "*?[") else [pat])

    known = set()
    for r in csv.DictReader(open(a.events), delimiter="\t"):
        known.add((r["insert_key"], r["context_key"]))
    print("existing events keyed: %d" % len(known))

    st = collections.Counter()
    seen_new = {}
    rows = []
    global _OPT
    _OPT = (a.min_insert, a.context, a.max_pairs_per_locus)
    if a.threads > 1:
        import multiprocessing
        with multiprocessing.Pool(a.threads) as pool:
            chunks = pool.map(_one_dir, dirs, chunksize=1)
    else:
        chunks = [_one_dir(d) for d in dirs]
    for cst, crows in chunks:
        st.update(cst)
        rows.extend(crows)

    # ---- STAGE 3, global: it dedups against the catalogue and against the
    # other new pairs, so it cannot be done per panel.
    emp, nst = [], []
    for r in rows:
        k = (r["insert_key"], r["context_key"])
        if k in known:
            st["already_in_catalogue"] += 1
            continue
        if k in seen_new:
            st["duplicate_of_another_new_pair"] += 1
            continue
        seen_new[k] = 1
        st["NOVEL_pairs"] += 1
        (emp if r["target_is_empty"] else nst).append(r)
    _write(a.out + "_empty_target_pairs.tsv", emp)
    _write(a.out + "_nested_target_pairs.tsv", nst)

    print()
    print("STAGE 1  loci with >=3 alleles          %8d" % st["loci_with_3plus_alleles"])
    print("STAGE 1  candidate non-extreme pairs    %8d" % st["candidate_pairs"])
    print("STAGE 2  decomposed at the SAME gates   %8d  (%.1f%% of candidates)"
          % (st["decomposed"], 100.0 * st["decomposed"] / max(1, st["candidate_pairs"])))
    print("         failed the gates               %8d" % st["failed_gates"])
    print("         context too short              %8d" % st["short_context"])
    print("STAGE 3  already an event in catalogue  %8d" % st["already_in_catalogue"])
    print("         duplicate of another new pair  %8d" % st["duplicate_of_another_new_pair"])
    print("STAGE 3  NOVEL pairs                    %8d" % st["NOVEL_pairs"])
    print("           empty target  (catalogue-comparable)  %8d" % len(emp))
    print("           nested target (NOT a bag flank)       %8d" % len(nst))
    print()
    print("wrote %s_empty_target_pairs.tsv and %s_nested_target_pairs.tsv"
          % (a.out, a.out))
    return 0


COLS = ["panel", "locus_id", "parent_locus", "short_allele", "long_allele",
        "short_len", "long_len", "inserted_len", "offset_in_short",
        "insert_md5", "insert_key", "context_key", "identity", "coverage",
        "target_lost", "overlap", "target_is_empty"]


def _write(path, rows):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w") as fh:
        fh.write("\t".join(COLS) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in COLS) + "\n")


def _one_dir(d):
    """Stages 1-2 for one panel. Returns (counter, rows). Pure: no globals
    written, so it is safe in a process pool."""
    min_insert, context, max_pairs = _OPT
    st, rows = collections.Counter(), []
    lp, ap_ = os.path.join(d, "loci.tsv"), os.path.join(d, "alleles.fna")
    if os.path.exists(lp) and os.path.exists(ap_):
        panel = LI.panel_label(d)
        aseq = LI.read_alleles(ap_)
        for r in csv.DictReader(open(lp), delimiter="\t"):
            av = aseq.get(r["locus_id"], {})
            if len(av) < 3:
                continue
            st["loci_with_3plus_alleles"] += 1
            short_a = r["shortest_allele"]
            long_a = r["longest_allele"]
            # dedupe IDENTICAL sequences before pairing: a locus with eight
            # alleles that are really three sequences must yield pairs from
            # the three, not from the eight.
            uniq = {}
            for aid, sq in av.items():
                uniq.setdefault(sq, aid)
            byl = sorted(uniq.items(), key=lambda kv: len(kv[0]))
            short_seq = byl[0][0] if byl else ""
            done = 0
            for i in range(len(byl)):
                if done >= max_pairs:
                    break
                for j in range(i + 1, len(byl)):
                    if done >= max_pairs:
                        break
                    s_seq, s_id = byl[i]
                    l_seq, l_id = byl[j]
                    if len(s_seq) == len(l_seq):
                        continue          # substitution only, cannot stack
                    if {s_id, l_id} == {short_a, long_a}:
                        continue          # 70_ already reported this one
                    st["candidate_pairs"] += 1
                    done += 1
                    td = L2.tolerant_decompose(s_seq, l_seq, 0.90, 0.90,
                                               min_insert, 3.0)
                    if td is None or not td["ok"]:
                        st["failed_gates"] += 1
                        continue
                    st["decomposed"] += 1
                    ins = td["insert_seq"]
                    off = td["lcp"]
                    lo = max(0, off - context)
                    hi = min(len(s_seq), off + context)
                    ctx = s_seq[lo:hi]
                    if len(ctx) < context:
                        st["short_context"] += 1
                        continue
                    rows.append({
                        "panel": panel, "locus_id": r["locus_id"],
                        "parent_locus": LI.qualify(panel, r["locus_id"]),
                        "short_allele": s_id, "long_allele": l_id,
                        "short_len": len(s_seq), "long_len": len(l_seq),
                        "inserted_len": len(ins), "offset_in_short": off,
                        "insert_md5": LI.md5(ins),
                        "insert_key": LI.canonical_insert_key(
                            l_seq, td["insert_start_long"], len(ins)),
                        "context_key": LI.canon_key(ctx),
                        "identity": "%.4f" % td["identity"],
                        "coverage": "%.4f" % td["coverage"],
                        "target_lost": td["target_lost"],
                        "overlap": td["ambiguity"],
                        # the target is the PRE-INSERTION site only when the
                        # shorter member of the pair is the locus's shortest
                        "target_is_empty": (s_seq == short_seq)})
    return st, rows


if __name__ == "__main__":
    sys.exit(main())
