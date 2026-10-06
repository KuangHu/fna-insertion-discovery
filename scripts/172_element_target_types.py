#!/usr/bin/env python3
"""What TYPE of mobile element, and what TYPE of gene does it land in?

171_ records, per event, which Pfam model fired on the element's own ORFs and
which fired on the disrupted host gene. This turns those into the two things
actually worth reading:

  1 element type   grouped into MGE classes from the Pfam family that fired,
                   so "IS3-family transposase" and "IS110 transposase" are
                   distinguishable rather than both being "transposase"
  2 target type    the Pfam family of the disrupted gene
  3 the cross-tab  which element classes land in which gene classes

MGE CLASSES, assigned from the Pfam family name. The groups are coarse on
purpose -- Pfam has dozens of transposase families and splitting them all
would give singleton rows. A family that matches nothing is reported under its
own name rather than being dropped into "other", so an unanticipated class
stays visible instead of disappearing into a bucket.

WHAT A TYPE IS. A Pfam hit above the gathering threshold says the protein
contains that domain. For the element that is strong evidence of class: a DDE
transposase domain is what makes an IS an IS. For the target it is weaker --
a domain names a function, not a gene, and most bacterial ORFs hit nothing at
all. `no_pfam_hit` is reported as its own row rather than hidden, because for
target genes it will be a large share and treating it as missing data would
bias every rate computed over the named ones.

SURVIVOR BIAS applies to the whole table. Insertions into essential genes
killed their host and were never sampled, so the target types seen are those
a cell can live without. Nothing here is a map of where elements insert; it
is a map of where elements insert AND the cell survives.
"""
import collections
import csv
import sys

csv.field_size_limit(sys.maxsize)

# element class <- substring of the Pfam family name, first match wins, so
# order matters: the specific families come before the generic ones.
MGE_CLASS = [
    ("IS110/IS492", ("is110", "transposase_20", "dedd_tnp")),
    ("IS200/IS605", ("is200", "is605", "transposase_17", "transposase_7")),
    ("IS3/IS150",   ("rve", "is3", "is150", "hth_transposase")),
    ("IS4/IS10",    ("is4", "is10", "transposase_11")),
    ("IS5/IS1182",  ("is5", "is1182", "transposase_mut")),
    ("IS21",        ("is21",)),
    ("IS30",        ("is30",)),
    ("IS66",        ("is66",)),
    ("IS91/helitron", ("is91", "helitron", "rep_3")),
    ("IS1",         ("is1",)),
    ("Tn3",         ("tn3", "tnpa")),
    ("Tn916/conj",  ("tn916", "relaxase", "mobc", "mob_pre", "traa")),
    ("serine recombinase", ("resolvase", "recombinase", "serine_recomb")),
    ("tyrosine integrase", ("phage_integr", "integrase", "int_", "xerc", "xerd")),
    ("DDE transposase", ("dde_", "transpos",)),
]


def mge_class(name):
    n = (name or "").lower()
    if not n or n == ".":
        return "no_pfam_hit"
    for cls, keys in MGE_CLASS:
        if any(k in n for k in keys):
            return cls
    return name          # keep it visible rather than bucketing it as "other"


rows = []
for p in sys.argv[1:]:
    rows.extend(csv.DictReader(open(p), delimiter="\t"))
if not rows:
    print("no rows")
    sys.exit(0)

mob = [r for r in rows if r.get("element_is_mobile") == "yes"]
print("intragenic events: %d   with an MGE-classified insert: %d (%.1f%%)"
      % (len(rows), len(mob), 100.0 * len(mob) / len(rows)))
print()

ec = collections.Counter(mge_class(r.get("element_pfam_name")) for r in mob)
print("ELEMENT TYPE")
print("  %-24s %7s" % ("class", "events"))
for k, v in ec.most_common():
    print("  %-24s %7d" % (k, v))
print()

tc = collections.Counter(r.get("target_pfam_name") or "." for r in mob)
named = sum(v for k, v in tc.items() if k not in (".", ""))
print("TARGET GENE TYPE  (%d of %d events have a Pfam hit, %.1f%%)"
      % (named, len(mob), 100.0 * named / max(1, len(mob))))
print("  %-40s %7s" % ("Pfam family of the disrupted gene", "events"))
for k, v in tc.most_common(25):
    print("  %-40s %7d" % (("no_pfam_hit" if k in (".", "") else k)[:40], v))
print()

print("CROSS-TAB  element class x target family (top pairs)")
print("  %-22s %-34s %6s" % ("element", "target", "n"))
x = collections.Counter((mge_class(r.get("element_pfam_name")),
                         r.get("target_pfam_name") or "no_pfam_hit")
                        for r in mob)
for (e, t), v in x.most_common(25):
    print("  %-22s %-34s %6d" % (e[:22], ("no_pfam_hit" if t == "." else t)[:34], v))
print()

print("Survivor bias: these are the genes a cell can live without having")
print("disrupted. Not a map of where elements insert -- a map of where they")
print("insert AND the host survives.")
