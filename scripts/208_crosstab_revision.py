#!/usr/bin/env python3
"""Cross-tabulate the revision and bound what may enter the official FASTA.

Three invariants are asserted, not merely printed, because each would be an
easy way to overstate the result:

  I1  official <= localization-resolved. Revising breakpoints must not promote
      a row whose LOCUS correspondence is still unresolved.
  I2  no row with candidate_search_truncated=yes is official, and none is
      called uniquely localized: a better copy may never have been examined.
  I3  every official row has a breakpoint_range, so a canonical coordinate is
      never presented as the unique biological cut point.

Counts are reported at the occurrence grain, with events and (event,carrier)
pairs given separately so the three are not conflated.
"""
import collections
import glob
import sys

SRC = sys.argv[1] if len(sys.argv) > 1 else \
    "/global/scratch/users/kh36969/fna_ins_discovery/inferred_empty_v2"


def rows(path):
    out = []
    files = sorted(glob.glob(path + "/parts/occ_*.tsv"))
    hdr = None
    import importlib.util, os
    sp = importlib.util.spec_from_file_location(
        "m206", os.path.dirname(os.path.abspath(__file__)) + "/206_revise_carrier_empty.py")
    m = importlib.util.module_from_spec(sp); sp.loader.exec_module(m)
    hdr = m.NEW_COLS
    for f in files:
        for line in open(f):
            x = line.rstrip("\n").split("\t")
            if len(x) == len(hdr):
                out.append(dict(zip(hdr, x)))
    return out


def shift_class(r):
    try:
        us, ds = int(r["up_shift"]), int(r["dn_shift"])
    except (ValueError, KeyError):
        return "no_shift_data"
    if us == 0 and ds == 0:
        return "no_shift"
    if us != ds:
        return "incoherent_sides"
    if r.get("shift_in_equiv_range") == "yes":
        return "within_equivalence_range"
    return "beyond_equivalence_range"


def main():
    R = rows(SRC)
    if not R:
        sys.exit("no rows under %s" % SRC)
    ev = set(r["db_id"] for r in R)
    pair = set((r["db_id"], r["carrier_accession"]) for r in R)
    print("=" * 72)
    print("GRAINS   events %d   (event,carrier) %d   occurrence rows %d"
          % (len(ev), len(pair), len(R)))

    ct = collections.Counter()
    for r in R:
        ct[(r["localization_class_revised"], shift_class(r),
            r["reconstructed_sequence_status"], r["official_fasta"])] += 1
    print("\nlocalization x shift_class x reconstructed_sequence x official")
    print("  %-26s %-24s %-12s %-8s %7s" % ("localization", "shift_class",
                                            "sequence", "official", "n"))
    for k, v in sorted(ct.items(), key=lambda x: -x[1]):
        print("  %-26s %-24s %-12s %-8s %7d" % (k[0], k[1], k[2], k[3], v))

    loc_res = sum(1 for r in R if r["localization_class_revised"] == "resolved")
    off = sum(1 for r in R if r["official_fasta"] == "yes")
    within = [r for r in R if shift_class(r) == "within_equivalence_range"]
    w_locres = sum(1 for r in within if r["localization_class_revised"] == "resolved")
    w_ver = sum(1 for r in within if r["sequence_equivalence_verified"] == "yes")
    added = sum(1 for r in R if "ADDED to official" in r.get("revision_note", ""))
    removed = sum(1 for r in R if "REMOVED from official" in r.get("revision_note", ""))
    trunc = [r for r in R if r.get("candidate_search_truncated") == "yes"]
    tr_off = sum(1 for r in trunc if r["official_fasta"] == "yes")
    tr_res = sum(1 for r in trunc if r["localization_class_revised"] == "resolved")
    bpu = collections.Counter(r["breakpoint_status"] for r in R)

    print("\nWITHIN-EQUIVALENCE-RANGE ROWS           %7d" % len(within))
    print("   of these, localization resolved       %7d" % w_locres)
    print("   of these, sequence equivalence VERIFIED %5d" % w_ver)
    print("\nOFFICIAL FASTA")
    print("   localization resolved (ceiling)       %7d" % loc_res)
    print("   official                              %7d" % off)
    print("   added vs v1                           %7d" % added)
    print("   removed vs v1                         %7d" % removed)
    print("\nBREAKPOINT STATUS")
    for k, v in bpu.most_common():
        print("   %-32s %7d" % (k, v))
    print("\nTRUNCATED SEARCHES                      %7d" % len(trunc))
    print("   called uniquely localized             %7d" % tr_res)
    print("   entered official FASTA                %7d" % tr_off)

    bad = []
    if off > loc_res:
        bad.append("I1 VIOLATED: official %d > localization-resolved %d" % (off, loc_res))
    if tr_off or tr_res:
        bad.append("I2 VIOLATED: truncated rows official=%d resolved=%d" % (tr_off, tr_res))
    miss = [r for r in R if r["official_fasta"] == "yes"
            and r.get("breakpoint_range_start", ".") in (".", "")]
    if miss:
        bad.append("I3 VIOLATED: %d official rows lack a breakpoint_range" % len(miss))
    print()
    if bad:
        for b in bad:
            print("  " + b)
        sys.exit(1)
    print("  I1 official <= localization-resolved      OK")
    print("  I2 no truncated row official or unique    OK")
    print("  I3 every official row has a range         OK")


if __name__ == "__main__":
    main()
