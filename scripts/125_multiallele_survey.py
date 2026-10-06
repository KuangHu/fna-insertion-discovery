#!/usr/bin/env python3
"""How much information is 70_ leaving on the table by using only the extremes?

70_ decomposes the LONGEST allele against the SHORTEST. Every other allele at
the locus is kept in alleles.tsv/alleles.fna and used only for stability
checks. Where a locus carries an intermediate allele, that pair

    A  L-------R          B  L--X----R          C  L--X--Y--R

is decomposable as A->B (X) and B->C (Y), whereas A->C reports one combined
X+Y difference. This survey does NOT decompose anything; it counts how often
the structure exists, so the value of building it can be judged before it is
built.

Reported per locus:
  n_alleles           distinct allele sequences at the anchor pair
  n_len_classes       distinct LENGTHS -- alleles of equal length differ by
                      substitution and cannot stack as nested insertions
  n_intermediate      length classes strictly between shortest and longest
  carriers_mid        genomes supporting the intermediate classes; an
                      intermediate seen in one genome is far weaker evidence
                      than one seen in many

An intermediate length class is necessary, not sufficient: whether A->B and
B->C actually decompose has to be tested. This bounds the opportunity above.
"""
import collections, csv, glob, os, sys
csv.field_size_limit(sys.maxsize)

tot = collections.Counter()
hist_all = collections.Counter()
hist_mid = collections.Counter()
mid_carrier = collections.Counter()
per_species = collections.defaultdict(lambda: [0, 0, 0])   # loci, with_mid, mid_classes

for d in sorted(glob.glob(sys.argv[1])):
    ap = os.path.join(d, "alleles.tsv")
    lp = os.path.join(d, "loci.tsv")
    if not (os.path.exists(ap) and os.path.exists(lp)):
        continue
    sp = os.path.basename(os.path.abspath(d).rstrip(os.sep).rsplit("/l2/", 1)[0])
    keep = set()
    for r in csv.DictReader(open(lp), delimiter="\t"):
        if r["event_class"] in ("insertion_target_retained", "replacement",
                                "length_polymorphism", "complex"):
            keep.add(r["locus_id"])
    by = collections.defaultdict(list)
    for r in csv.DictReader(open(ap), delimiter="\t"):
        if r["locus_id"] in keep:
            by[r["locus_id"]].append((int(r["length"]), int(r["n_genomes"])))
    for lid, al in by.items():
        lens = sorted({L for L, _ in al})
        n_mid = max(0, len(lens) - 2)
        cm = sum(g for L, g in al if lens and lens[0] < L < lens[-1])
        tot["loci"] += 1
        hist_all[min(len(al), 10)] += 1
        hist_mid[min(n_mid, 10)] += 1
        per_species[sp][0] += 1
        if n_mid:
            tot["loci_with_intermediate"] += 1
            tot["intermediate_classes"] += n_mid
            per_species[sp][1] += 1
            per_species[sp][2] += n_mid
            mid_carrier[min(cm, 10)] += 1

n = max(1, tot["loci"])
print("loci surveyed                     %8d" % tot["loci"])
print("loci with >=1 INTERMEDIATE length %8d  (%.1f%%)"
      % (tot["loci_with_intermediate"], 100.0 * tot["loci_with_intermediate"] / n))
print("total intermediate length classes %8d" % tot["intermediate_classes"])
print()
print("alleles per locus:")
for k in sorted(hist_all):
    print("  %2d%s : %7d" % (k, "+" if k == 10 else " ", hist_all[k]))
print()
print("intermediate length classes per locus:")
for k in sorted(hist_mid):
    print("  %2d%s : %7d" % (k, "+" if k == 10 else " ", hist_mid[k]))
print()
print("genomes supporting the intermediate class(es):")
for k in sorted(mid_carrier):
    print("  %2d%s genomes : %7d loci" % (k, "+" if k == 10 else " ", mid_carrier[k]))
print()
print("%-22s %9s %9s %7s %9s" % ("run", "loci", "with_mid", "pct", "mid_cls"))
for sp in sorted(per_species, key=lambda x: -per_species[x][1]):
    a, b, c = per_species[sp]
    print("%-22s %9d %9d %6.1f%% %9d" % (sp, a, b, 100.0 * b / max(1, a), c))
