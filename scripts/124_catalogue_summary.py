#!/usr/bin/env python3
"""Per-species summary of a catalogue directory: the numbers result_note quotes."""
import csv, glob, os, sys, collections
csv.field_size_limit(sys.maxsize)

print("%-16s %8s %8s %8s %7s %8s %8s %9s %8s" %
      ("species", "loci", "EVENTS", "redund", "S1%", "mob%", "both%",
       "frame_rep", "unres"))
tl = te = ts = tm = tb = tu = 0
tboth = 0
for p in sorted(glob.glob(sys.argv[1] + "/*_loci.tsv")):
    sp = os.path.basename(p).replace("_loci.tsv", "")
    n = s1 = mo = bad = unres = 0
    evs, ev_both = set(), set()
    for r in csv.DictReader(open(p), delimiter="\t"):
        n += 1
        evs.add(r["event_id"])
        a = r["S1_structurally_clean"] == "True"
        b = r.get("mobility_positive") == "True"
        s1 += a; mo += b
        if a and b:
            ev_both.add(r["event_id"])
        f = r.get("insert_frame_status", "")
        bad += f == "repaired"
        unres += f == "unresolved" or not r.get("insert_seq")
    ne = len(evs - {"."})
    nb = len(ev_both - {"."})
    print("%-16s %8d %8d %7.2fx %6.1f%% %7.1f%% %7.1f%% %9d %8d"
          % (sp, n, ne, n / max(1, ne), 100.0*s1/max(1,n), 100.0*mo/max(1,n),
             100.0*nb/max(1,ne), bad, unres))
    tl += n; te += ne; ts += s1; tm += mo; tb += bad; tu += unres; tboth += nb
print("%-16s %8d %8d %7.2fx %6.1f%% %7.1f%% %7.1f%% %9d %8d"
      % ("TOTAL", tl, te, tl/max(1,te), 100.0*ts/max(1,tl), 100.0*tm/max(1,tl),
         100.0*tboth/max(1,te), tb, tu))
