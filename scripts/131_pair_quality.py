#!/usr/bin/env python3
"""Do the multi-allele pairs look like the catalogue, or like recombination?

130_ found that the new-site yield per species tracks something other than
insertion biology: H. pylori returns 4,785 new sites against 1,645 existing
events (291%) and N. gonorrhoeae 181%, while E. coli returns 28.5%. Both
outliers are naturally competent, heavily recombinogenic species whose loci
carry many divergent haplotypes of differing length.

If an "intermediate allele" in those species is a recombinant haplotype rather
than a nested insertion state, then decomposing it against the shortest allele
manufactures an "insert" out of a divergent segment. That decomposition would
pass the gates, but it should sit CLOSER TO THEM than a real element does.

So compare the pairs to the gate they had to clear (identity >= 0.90):

  id_p10 / id_med   a clean nested insertion is near-identical outside the
                    inserted block, so identity should sit at ~1.00 and the
                    10th percentile well above the gate
  at_gate%          share below 0.93, i.e. within 3 points of the threshold
  cov_med           coverage of the short allele by the alignment

This is a distributional comparison ACROSS species, not a threshold. It says
which species' yield is suspect; it does not by itself label any single pair
wrong.
"""
import csv
import glob
import os
import sys

csv.field_size_limit(sys.maxsize)

W = sys.argv[1]
print("%-16s %8s %8s %8s %8s %9s %9s"
      % ("species", "n", "id_p10", "id_med", "cov_med", "len_med", "at_gate%"))
rows = []
for p in sorted(glob.glob(os.path.join(W, "multiallele",
                                       "*_empty_target_pairs.tsv"))):
    sp = os.path.basename(p).replace("_empty_target_pairs.tsv", "")
    if sp.startswith("v2_"):
        continue
    ide, cov, ln = [], [], []
    for r in csv.DictReader(open(p), delimiter="\t"):
        ide.append(float(r["identity"]))
        cov.append(float(r["coverage"]))
        ln.append(int(r["inserted_len"]))
    if not ide:
        continue
    ide.sort()
    cov.sort()
    ln.sort()
    gate = 100.0 * sum(1 for x in ide if x < 0.93) / len(ide)
    rows.append((sp, len(ide), ide[len(ide) // 10], ide[len(ide) // 2],
                 cov[len(cov) // 2], ln[len(ln) // 2], gate))
for r in sorted(rows, key=lambda x: -x[6]):
    print("%-16s %8d %8.4f %8.4f %8.4f %9d %8.1f%%" % r)
print()
print("Sorted by at_gate%. A species whose pairs cluster near the identity")
print("gate is reporting divergence it had to strain to call an insertion.")
