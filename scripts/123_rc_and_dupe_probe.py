#!/usr/bin/env python3
"""Finding #2, measured two ways, without the sampling bias of 122_.

122_ ran its revcomp test only on loci whose re-slice was already wrong -- a
subset selected for having upstream indels, i.e. for messy alignments. 87.8%
disagreement there is not a rate for the corpus.

TEST 1 (mechanism, unbiased): random decomposable loci of BOTH methods, each
  re-decomposed in the forward and reverse-complement frame. Does the junction
  the algorithm picks survive reverse-complementation?

TEST 2 (consequence, the number that matters): does the mechanism actually
  split real events in the shipped dedup output? For every event, key its
  insert two ways --
      exact    canon(insert)                     <- what 100_ uses
      trimmed  canon(insert[T:-T])               <- immune to junction jitter
                                                     of up to T bp at each end
  and, holding the context key fixed, count events that merge under `trimmed`
  but were emitted as separate event_ids. That is double-counting, measured on
  the output rather than argued from a synthetic case.
"""
import collections, csv, glob, hashlib, importlib.util, os, random, sys
csv.field_size_limit(sys.maxsize)
D = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, D)
import lib_insert as LI
spec = importlib.util.spec_from_file_location(
    "l2", os.path.join(D, "70_allele_reconstructor.py"))
l2 = importlib.util.module_from_spec(spec); spec.loader.exec_module(l2)

TRIM = 30
rng = random.Random(23)

rows = []
for d in sorted(glob.glob(sys.argv[1])):
    lp, ap = os.path.join(d, "loci.tsv"), os.path.join(d, "alleles.fna")
    if not (os.path.exists(lp) and os.path.exists(ap)):
        continue
    aseq = LI.read_alleles(ap)
    for r in csv.DictReader(open(lp), delimiter="\t"):
        if r["event_class"] not in ("insertion_target_retained", "replacement"):
            continue
        av = aseq.get(r["locus_id"], {})
        ls, ss = av.get(r["longest_allele"], ""), av.get(r["shortest_allele"], "")
        if not ls or not ss:
            continue
        ins, frame = LI.extract_insert(r, ls)
        if not ins:
            continue
        rows.append({"lid": r["locus_id"], "ins": ins, "short": ss, "long": ls,
                     "lcp": int(r["lcp_bp"]), "method": r["decomposition_method"],
                     "overlap": int(r.get("junction_ambiguity_bp") or 0),
                     "frame": frame})
print("decomposable loci with a verified insert: %d" % len(rows))

# ---------------- TEST 1 ----------------
print("\n=== TEST 1: is the junction choice revcomp-invariant? (unbiased) ===")
by = collections.defaultdict(lambda: [0, 0])      # method -> [n, invariant]
ov = collections.defaultdict(lambda: [0, 0])      # overlap band -> [n, invariant]
sample = rows if len(rows) <= 3000 else rng.sample(rows, 3000)
for r in sample:
    ss, ls = r["short"], r["long"]
    f = l2.tolerant_decompose(ss, ls, 0.90, 0.90, 500, 3.0)
    v = l2.tolerant_decompose(LI.rc(ss), LI.rc(ls), 0.90, 0.90, 500, 3.0)
    if not f or not v or not f["ok"] or not v["ok"]:
        continue
    same = LI.canon_key(f["insert_seq"]) == LI.canon_key(v["insert_seq"])
    by[r["method"]][0] += 1; by[r["method"]][1] += same
    b = ("0" if r["overlap"] == 0 else "1-3" if r["overlap"] <= 3
         else "4-15" if r["overlap"] <= 15 else ">15")
    ov[b][0] += 1; ov[b][1] += same
print("%-22s %8s %10s %8s" % ("recorded method", "n", "invariant", "pct"))
for k in sorted(by):
    n, s = by[k]; print("%-22s %8d %10d %7.1f%%" % (k, n, s, 100.0*s/max(1,n)))
n = sum(v[0] for v in by.values()); s = sum(v[1] for v in by.values())
print("%-22s %8d %10d %7.1f%%" % ("ALL", n, s, 100.0*s/max(1,n)))
print("\nby junction overlap (the direct-repeat axis the finding names):")
for k in ("0", "1-3", "4-15", ">15"):
    if k in ov:
        a, b_ = ov[k]; print("  overlap %-5s %7d %9d %7.1f%%" % (k, a, b_, 100.0*b_/max(1,a)))

# ---------------- TEST 2 ----------------
print("\n=== TEST 2: does it split real events in the shipped output? ===")
ev_of = {}
for r in csv.DictReader(open(sys.argv[2]), delimiter="\t"):
    ev_of[r["locus_id"]] = r["event_id"]
ctx_of = {}
for r in csv.DictReader(open(sys.argv[2]), delimiter="\t"):
    ctx_of[r["locus_id"]] = r["context_key"]

# one representative insert per event
rep = {}
for r in rows:
    e = ev_of.get(r["lid"])
    if e and e not in rep:
        rep[e] = (r["ins"], ctx_of.get(r["lid"], ""))
print("events with a representative insert: %d" % len(rep))

exact_g = collections.defaultdict(set)
trim_g = collections.defaultdict(set)
for e, (ins, ctx) in rep.items():
    exact_g[(ctx, LI.canon_key(ins))].add(e)
    t = ins[TRIM:-TRIM] if len(ins) > 2 * TRIM + 20 else ins
    trim_g[(ctx, LI.canon_key(t))].add(e)
split_exact = sum(1 for g in exact_g.values() if len(g) > 1)
merged = [g for g in trim_g.values() if len(g) > 1]
extra = sum(len(g) - 1 for g in merged)
print("groups of >1 event sharing context+EXACT insert key : %d" % split_exact)
print("groups of >1 event sharing context+TRIMMED(%dbp) key: %d" % (TRIM, len(merged)))
print("events that would collapse under the trimmed key    : %d of %d (%.2f%%)"
      % (extra, len(rep), 100.0 * extra / max(1, len(rep))))
for g in sorted(merged, key=lambda x: -len(x))[:5]:
    print("   would merge: %s" % ", ".join(sorted(g)))
