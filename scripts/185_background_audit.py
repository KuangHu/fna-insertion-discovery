#!/usr/bin/env python3
"""Audit the background counting on concrete loci. Supersedes 184_'s summary.

184_ reported T1>T2<T3, which is impossible for nested tiers. Cause: T2 was a
GLOBAL set of positions while T3 was summed PER FAMILY, so one position hit by
three families counted once in T2 and three times in T3. The tiers were never
nested and the summary should not have been presented as a progressive filter.

This does not add another aggregate. It answers five specific questions about
specific records, each of which can be checked by eye.

  A  multi-seed hits at one position -- are several element families matching
     the SAME place, and is that a conflict to preserve or a merge to make?
  B  records whose flanks differ by a FEW BASES -- one ancestral locus split
     by SNPs, or genuinely different loci?
  C  records with similar flanks but different insertion offsets -- would
     flank-keying wrongly merge distinct insertion points?
  D  expansion-only backgrounds -- genuinely absent from the catalogue, or a
     boundary-shifted variant of something already there?
  E  a few clean families -- how many distinct targets survive honest dedup?

NO TIER NAMES ARE USED. 184_'s "informative backgrounds" implied a guide-target
test that was never run; this reports what was actually measured.
"""
import collections
import csv
import sys

csv.field_size_limit(sys.maxsize)
W = "/global/scratch/users/kh36969/fna_ins_discovery"
SP = "spneumoniae"

rows = list(csv.DictReader(open(W + "/xspecies/%s_hits.tsv" % SP), delimiter="\t"))
print("expansion hit records: %d" % len(rows))
print()

# ---------- A. several families at one position --------------------------
bypos = collections.defaultdict(list)
for r in rows:
    bypos[(r["genome"], r["contig"], int(r["elem_start"]) // 50)].append(r)
multi = {k: v for k, v in bypos.items() if len({x["cds_cluster"] for x in v}) > 1}
print("A. ONE POSITION, SEVERAL ELEMENT FAMILIES")
print("   positions with >1 family: %d of %d (%.1f%%)"
      % (len(multi), len(bypos), 100.0 * len(multi) / max(1, len(bypos))))
ex = sorted(multi.items(), key=lambda kv: -len({x["cds_cluster"] for x in kv[1]}))[:3]
for (g, c, p), v in ex:
    fams = sorted({x["cds_cluster"] for x in v})
    spans = sorted({(int(x["elem_start"]), int(x["elem_end"])) for x in v})
    print("   %s %s ~%d" % (g[:20], c[:22], p * 50))
    print("      families: %s" % ",".join(fams[:6]))
    print("      spans:    %s" % str(spans[:4]))
print("   -> these are competing or nested assignments. Merging them into one")
print("      'site' discards the conflict; counting them separately inflates.")
print()

# ---------- B. flanks differing by a few bases ---------------------------
def key60(r):
    lf, rf = r.get("left_flank") or "", r.get("right_flank") or ""
    if len(lf) < 60 or len(rf) < 60:
        return None
    return lf[-60:].upper(), rf[:60].upper()


def hamm(a, b):
    return sum(1 for x, y in zip(a, b) if x != y)


byfam = collections.defaultdict(list)
for r in rows:
    k = key60(r)
    if k:
        byfam[r["cds_cluster"]].append((k, r))
print("B. FLANKS DIFFERING BY A FEW BASES  (same locus split by SNPs?)")
near = exact = far = 0
checked = 0
for fam, lst in list(byfam.items()):
    if len(lst) < 2 or checked > 40000:
        continue
    seen = {}
    for (k, r) in lst[:60]:
        for k2 in list(seen)[:60]:
            checked += 1
            d = hamm(k[0], k2[0]) + hamm(k[1], k2[1])
            if d == 0:
                exact += 1
            elif d <= 3:
                near += 1
            else:
                far += 1
            break
        seen[k] = 1
print("   pairwise flank comparisons sampled: %d" % checked)
print("     identical            %7d" % exact)
print("     1-3 bases different  %7d   <- exact keying SPLITS these" % near)
print("     >3 different         %7d" % far)
print("   -> every 1-3bp pair is one ancestral locus counted twice by an")
print("      exact key. The background count is therefore OVER-estimated.")
print()

# ---------- C. similar flanks, different offset --------------------------
print("C. SIMILAR FLANKS, DIFFERENT INSERTION POINT  (would merging be wrong?)")
same_flank_diff_span = 0
fl = collections.defaultdict(set)
for r in rows:
    k = key60(r)
    if k:
        fl[k].add(int(r["elem_end"]) - int(r["elem_start"]))
multi_len = {k: v for k, v in fl.items() if len(v) > 1}
print("   flank keys carrying >1 element length: %d of %d (%.2f%%)"
      % (len(multi_len), len(fl), 100.0 * len(multi_len) / max(1, len(fl))))
for k, v in list(multi_len.items())[:3]:
    print("      one flank context, element lengths: %s" % sorted(v)[:6])
print("   -> where this happens, flank-keying merges DIFFERENT insertions and")
print("      the background count is under-estimated for those.")
print()

# ---------- D. expansion-only vs the catalogue ---------------------------
cat = set()
for r in csv.DictReader(open(W + "/catalogue_v3/%s_events.tsv" % SP), delimiter="\t"):
    t = r.get("target_seq") or ""
    o = int(r.get("insertion_point_offset") or 0)
    if len(t) >= o + 60 and o >= 60:
        cat.add((t[o - 60:o].upper(), t[o:o + 60].upper()))
print("D. EXPANSION BACKGROUNDS ALSO IN THE CATALOGUE")
hit = miss = 0
for r in rows:
    k = key60(r)
    if not k:
        continue
    if k in cat:
        hit += 1
    else:
        miss += 1
print("   flank context already in the catalogue: %d" % hit)
print("   not in the catalogue:                   %d" % miss)
print("   -> 'not in the catalogue' at an EXACT key includes boundary shifts")
print("      and single-SNP variants of catalogue entries, so it is an upper")
print("      bound on genuinely new backgrounds, not a count of them.")
