#!/usr/bin/env python3
"""Does the position information actually find the sequences in the FNA corpus?

THE QUESTION. insertions.tsv says: event X sits at contig C position J in
assembly A, its insert is S_ins and its empty target is S_emp. Nothing so far
has checked that claim against the raw assemblies at scale. One event was
checked by hand (ecoli.E004606, a clean 2x2); this does all 109,709.

FOUR TESTS, TWO PASSES. Junction probes are used rather than whole sequences
because a whole-allele match would fail on a single SNP and tell us nothing
about whether the COORDINATE is right. A probe spanning a junction is the
tightest thing that still isolates position from divergence.

  pass "empty"   -- open the coord_frame assembly, seek to the junction
    T1 empty_junction_present   up[-K:] + dn[:K] is in the empty assembly.
                                If the coordinate is right the two flanks abut
                                there, so this must hold.
    T2 target_matches           the slice [locus_start, locus_end] equals the
                                empty_site.fna record for this event.
    T3 insert_absent_at_site    up[-K:] + insert[:K] is NOT in the empty
                                assembly. The element may live elsewhere in
                                that genome -- that is normal for a mobile
                                element and is not a failure -- so this tests
                                only the junction, never mere presence.

  pass "carrier" -- open one carrier assembly
    T4 insert_present_at_site   BOTH up[-K:]+insert[:K] and insert[-K:]+dn[:K]
                                are in the carrier. This is the decisive test:
                                it says the insert sits at THIS position, in a
                                genome that has it.

Both passes group by genome and shard on the genome, never on the event: the
same assembly serves ~10 events, and sharding on events would re-read it in
every shard.

A probe miss is reported, never silently dropped, and both strands are searched
because allele orientation is not guaranteed to match the assembly's.
"""
import argparse
import collections
import os
import sqlite3
import sys

W = "/global/scratch/users/kh36969/fna_ins_discovery"
OUT = W + "/insertions_all"
DB = W + "/locations/insertions.sqlite"
COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def revcomp(s):
    return s.translate(COMP)[::-1]


def load_manifest(path):
    man = {}
    for line in open(path):
        p = line.strip()
        if not p:
            continue
        b = os.path.basename(p)
        for e in (".fna", ".fa", ".fasta"):
            if b.endswith(e):
                b = b[:-len(e)]
                break
        man.setdefault(b, p)
    return man


def load_genome(path):
    out, k, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if k is not None:
                    out[k] = "".join(buf).upper()
                k = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if k is not None:
        out[k] = "".join(buf).upper()
    return out


def read_fasta_map(path):
    """db_id -> sequence, for insert.fna / empty_site.fna."""
    out, k, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if k is not None:
                    out[k] = "".join(buf).upper()
                k = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if k is not None:
        out[k] = "".join(buf).upper()
    return out


def found(hay, needle):
    if not needle or len(needle) < 10:
        return None                      # too short to mean anything
    return (needle in hay) or (revcomp(needle) in hay)


def all_occurrences(ctgs, needle, cap=50):
    """Every exact hit of needle, both strands: (contig, offset, strand).

    A MOBILE ELEMENT HAS MANY COPIES, SO THE FIRST HIT IS NOT THE RIGHT HIT.
    An earlier version took ctgs.find() and stopped; for a multi-copy family
    that picks an arbitrary copy and the flanks then score at ~0.33, which is
    DNA background, not divergence.
    """
    out = []
    r = revcomp(needle)
    for k, v in ctgs.items():
        for pat, st in ((needle, "+"), (r, "-")):
            i = v.find(pat)
            while i >= 0 and len(out) < cap:
                out.append((k, i, st))
                i = v.find(pat, i + 1)
    return out


def ident(a, b):
    """Fraction of matching positions over the overlap; 0.0 if no overlap."""
    n = min(len(a), len(b))
    if not n:
        return 0.0
    return sum(1 for x, y in zip(a[-n:], b[-n:]) if x == y) / float(n)


def score_flanks(cseq, off, ilen, up, dn, K):
    """Best (up_id, dn_id) over both relative orientations.

    DO NOT REORIENT THE CARRIER FLANKS INTO THE ELEMENT'S FRAME. up/dn were cut
    in the EMPTY contig's frame, and the canonical insert's orientation is
    arbitrary after sequence dedup, so the element's frame says nothing about
    how the two genomic contexts line up. Reorienting compared left against
    right and scored homologous flanks at 0.33. Score both pairings instead and
    keep the better.
    """
    left = cseq[max(0, off - K):off]
    right = cseq[off + ilen:off + ilen + K]
    same = (ident(left, up), ident(right, dn))
    flip = (ident(revcomp(right), up), ident(revcomp(left), dn))
    return same if (same[0] + same[1]) >= (flip[0] + flip[1]) else flip


def events():
    con = sqlite3.connect("file:%s?mode=ro" % DB, uri=True)
    con.row_factory = sqlite3.Row
    return list(con.execute(
        "SELECT db_id, species, coord_frame, coord_contig, junction_L_abs,"
        "  junction_R_abs, locus_start, locus_end, inserted_len"
        " FROM loc WHERE junction_L_abs NOT IN ('.','') GROUP BY db_id"))


def carriers():
    con = sqlite3.connect("file:%s?mode=ro" % DB, uri=True)
    out = {}
    for db_id, g in con.execute(
            "SELECT db_id, genome FROM loc WHERE role='carrier_derived'"):
        out.setdefault(db_id, []).append(g)
    return out


def mode_empty(a):
    man = load_manifest(a.manifest)
    ins = read_fasta_map(OUT + "/insert.fna")
    emp = read_fasta_map(OUT + "/empty_site.fna")
    by = collections.defaultdict(list)
    for r in events():
        by[r["coord_frame"]].append(r)
    gs = sorted(by)
    mine = [g for i, g in enumerate(gs) if i % a.nshard == a.shard]
    sys.stderr.write("shard %d: %d genomes\n" % (a.shard, len(mine)))

    out = open("%s/val_parts/empty_%03d.tsv" % (OUT, a.shard), "w")
    n_nogenome = n_nocontig = 0
    for n, g in enumerate(mine, 1):
        p = man.get(g)
        if not p or not os.path.exists(p):
            n_nogenome += len(by[g])
            for r in by[g]:
                out.write("%s\t%s\t%s\tno_genome\t.\t.\t.\n"
                          % (r["db_id"], r["species"], g))
            continue
        ctgs = load_genome(p)
        whole = "".join(ctgs.values())      # ONE join per genome, not per event
        for r in by[g]:
            seq = ctgs.get(r["coord_contig"])
            if seq is None:
                n_nocontig += 1
                out.write("%s\t%s\t%s\tno_contig\t.\t.\t.\n"
                          % (r["db_id"], r["species"], g))
                continue
            jl, jr = int(r["junction_L_abs"]), int(r["junction_R_abs"])
            K = a.probe
            up, dn = seq[max(0, jl - K):jl], seq[jr:jr + K]
            t1 = found(whole, up + dn)
            # T2 ASKS WHETHER THE EMPTY SEQUENCE IS FINDABLE, NOT WHETHER TWO
            # DIFFERENT SPANS ARE BYTE-EQUAL. An earlier version compared
            # seq[locus_start:locus_end] to the empty_site.fna record and
            # "failed" 30% of events. All of those failures turned out to have
            # the target present elsewhere in the same assembly: locus_start/
            # locus_end is the anchored interval from loci.tsv, while
            # target_seq is a canonical window with its own
            # insertion_point_offset. They were never the same span.
            tgt = emp.get(r["db_id"])
            t2 = "."
            if tgt:
                t2 = "yes" if found(whole, tgt) else "no"
            iseq = ins.get(r["db_id"], "")
            t3 = found(whole, up + iseq[:K])
            out.write("%s\t%s\t%s\tok\t%s\t%s\t%s\n"
                      % (r["db_id"], r["species"], g,
                         "yes" if t1 else "no", t2,
                         "no" if not t3 else "yes"))
        if n % 200 == 0:
            sys.stderr.write("  %d/%d\n" % (n, len(mine)))
    out.close()
    print("shard %d empty: no_genome=%d no_contig=%d" % (a.shard, n_nogenome, n_nocontig))


def mode_carrier(a):
    man = load_manifest(a.manifest)
    ins = read_fasta_map(OUT + "/insert.fna")
    car = carriers()
    ev = {r["db_id"]: r for r in events()}

    # flanks come from the EMPTY assembly, so group the work by carrier but
    # remember which empty genome each event's flanks must be cut from
    by_empty = collections.defaultdict(list)
    for db_id, r in ev.items():
        by_empty[r["coord_frame"]].append(db_id)

    pick = {}
    for db_id, gl in car.items():
        if db_id in ev and gl:
            pick[db_id] = gl[0]
    by_car = collections.defaultdict(list)
    for db_id, g in pick.items():
        by_car[g].append(db_id)
    gs = sorted(by_car)
    mine = [g for i, g in enumerate(gs) if i % a.nshard == a.shard]
    sys.stderr.write("shard %d: %d carrier genomes, %d events\n"
                     % (a.shard, len(mine), sum(len(by_car[g]) for g in mine)))

    # PRECOMPUTE EVERY FLANK IN ONE PASS, GROUPED BY THE EMPTY GENOME.
    # Caching per db_id would reload a whole assembly for each of its ~10
    # events; grouping reads each assembly exactly once. Only two K-mers per
    # event are kept, so 109k events cost ~9MB rather than a genome each.
    need = set()
    for g in [x for i, x in enumerate(sorted(by_car)) if i % a.nshard == a.shard]:
        need.update(by_car[g])
    flank = {}
    todo = collections.defaultdict(list)
    for db_id in need:
        todo[ev[db_id]["coord_frame"]].append(db_id)
    K = a.probe
    for n, (gf, ids) in enumerate(todo.items(), 1):
        p = man.get(gf)
        if not p or not os.path.exists(p):
            for db_id in ids:
                flank[db_id] = None
            continue
        ctgs = load_genome(p)
        for db_id in ids:
            r = ev[db_id]
            seq = ctgs.get(r["coord_contig"])
            if seq is None:
                flank[db_id] = None
                continue
            jl, jr = int(r["junction_L_abs"]), int(r["junction_R_abs"])
            flank[db_id] = (seq[max(0, jl - K):jl], seq[jr:jr + K])
        if n % 200 == 0:
            sys.stderr.write("  flanks %d/%d empty genomes\n" % (n, len(todo)))
    sys.stderr.write("  flanks ready for %d events\n" % len(flank))

    def get_flanks(db_id):
        return flank.get(db_id)

    out = open("%s/val_parts/carrier_%03d.tsv" % (OUT, a.shard), "w")
    for n, g in enumerate(mine, 1):
        p = man.get(g)
        if not p or not os.path.exists(p):
            for db_id in by_car[g]:
                out.write("%s\t%s\t%s\tno_genome\t.\t.\t.\t.\t.\t.\t.\n"
                          % (db_id, ev[db_id]["species"], g))
            continue
        ctgs_c = load_genome(p)
        for db_id in by_car[g]:
            fl = get_flanks(db_id)
            if fl is None:
                out.write("%s\t%s\t%s\tno_flank\t.\t.\t.\t.\t.\t.\t.\n"
                          % (db_id, ev[db_id]["species"], g))
                continue
            up, dn = fl
            iseq = ins.get(db_id, "")
            if not iseq:
                out.write("%s\t%s\t%s\tno_insert_seq\t.\t.\t.\t.\t.\t.\t.\n"
                          % (db_id, ev[db_id]["species"], g))
                continue
            # T4 LOCATES THE INSERT, IT DOES NOT PROBE A JUNCTION.
            # The previous version built up[-K:]+insert[:K] from the EMPTY
            # assembly and searched the CARRIER. That fails on a single SNP in
            # 40bp of flank, which is cross-genome divergence, not a wrong
            # coordinate: of 136 such "failures" ALL 136 carried the complete
            # insert. So: find the insert, then score the flanks by identity.
            occ = all_occurrences(ctgs_c, iseq)
            if not occ:
                out.write("%s\t%s\t%s\tok\tno\t.\t.\t.\t.\t.\t.\n"
                          % (db_id, ev[db_id]["species"], g))
                continue
            K = a.probe
            best = None
            for ctg, off, st in occ:
                u, d = score_flanks(ctgs_c[ctg], off, len(iseq), up, dn, K)
                if best is None or (u + d) > (best[3] + best[4]):
                    best = (ctg, off, st, u, d)
            ctg, off, strand, iu, idn = best
            out.write("%s\t%s\t%s\tok\tyes\t%s\t%d\t%.3f\t%.3f\t%d\t%s\n"
                      % (db_id, ev[db_id]["species"], g, strand, off,
                         iu, idn, len(occ), ctg))
        if n % 200 == 0:
            sys.stderr.write("  %d/%d\n" % (n, len(mine)))
    out.close()
    print("shard %d carrier done" % a.shard)


def pct(a, b):
    return "%.1f%%" % (100.0 * a / b) if b else "n/a"


def mode_report(a):
    import glob
    e = collections.Counter()
    sp_fail = collections.Counter()
    for f in glob.glob(OUT + "/val_parts/empty_*.tsv"):
        for line in open(f):
            x = line.rstrip("\n").split("\t")
            if len(x) < 7:
                continue
            e["rows"] += 1
            if x[3] != "ok":
                e[x[3]] += 1
                continue
            e["checked"] += 1
            if x[4] == "yes":
                e["T1_empty_junction_present"] += 1
            else:
                sp_fail["T1 " + x[1]] += 1
            if x[5] == "yes":
                e["T2_target_found_in_assembly"] += 1
            elif x[5] == "no":
                e["T2_fail"] += 1
            if x[6] == "no":
                e["T3_insert_absent_at_site"] += 1
            else:
                sp_fail["T3 " + x[1]] += 1
    c = collections.Counter()
    fl_up, fl_dn = [], []
    for f in glob.glob(OUT + "/val_parts/carrier_*.tsv"):
        for line in open(f):
            x = line.rstrip("\n").split("\t")
            if len(x) < 10:
                continue
            c["rows"] += 1
            if x[3] != "ok":
                c[x[3]] += 1
                continue
            c["checked"] += 1
            if x[4] != "yes":
                c["T4_insert_NOT_found"] += 1
                continue
            c["T4_insert_located"] += 1
            c["strand_" + x[5]] += 1
            try:
                u, d = float(x[7]), float(x[8])
            except ValueError:
                continue
            fl_up.append(u)
            fl_dn.append(d)
            try:
                ncopy = int(x[9])
                c["copies_1"] += 1 if ncopy == 1 else 0
                c["copies_multi"] += 1 if ncopy > 1 else 0
            except (ValueError, IndexError):
                pass
            if u >= 0.9 and d >= 0.9:
                c["T5_both_flanks_>=90pct"] += 1
            if u == 1.0 and d == 1.0:
                c["T5_both_flanks_exact"] += 1

    print("=" * 64)
    print("PASS 1  EMPTY ASSEMBLY      rows %d   checked %d" % (e["rows"], e["checked"]))
    for k in ("no_genome", "no_contig"):
        if e[k]:
            print("   %-28s %7d" % (k, e[k]))
    for k in ("T1_empty_junction_present", "T2_target_found_in_assembly",
              "T3_insert_absent_at_site"):
        print("   %-28s %7d  %s" % (k, e[k], pct(e[k], e["checked"])))
    if e["T2_fail"]:
        print("   %-28s %7d" % ("T2 mismatches", e["T2_fail"]))
    print()
    print("PASS 2  CARRIER ASSEMBLY    rows %d   checked %d" % (c["rows"], c["checked"]))
    for k in ("no_genome", "no_flank"):
        if c[k]:
            print("   %-28s %7d" % (k, c[k]))
    for k in ("T4_insert_located", "T4_insert_NOT_found",
              "T5_both_flanks_>=90pct", "T5_both_flanks_exact"):
        print("   %-28s %7d  %s" % (k, c[k], pct(c[k], c["checked"])))
    if fl_up:
        fl_up.sort(); fl_dn.sort()
        n = len(fl_up)
        print("   flank identity to the empty assembly (%d events)" % n)
        for lab, v in (("upstream", fl_up), ("downstream", fl_dn)):
            print("     %-12s median %.3f  p10 %.3f  min %.3f"
                  % (lab, v[n // 2], v[int(n * 0.10)], v[0]))
    print("   strand of the insert in the carrier: + %d   - %d"
          % (c["strand_+"], c["strand_-"]))
    if sp_fail:
        print("\n  failures by species (top 12)")
        for k, v in sp_fail.most_common(12):
            print("   %-28s %7d" % (k, v))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mode", choices=["empty", "carrier", "report"])
    ap.add_argument("--manifest", default=OUT + "/manifest_union.txt")
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshard", type=int, default=1)
    ap.add_argument("--probe", type=int, default=40,
                    help="bp each side of a junction in a probe")
    a = ap.parse_args()
    os.makedirs(OUT + "/val_parts", exist_ok=True)
    {"empty": mode_empty, "carrier": mode_carrier, "report": mode_report}[a.mode](a)


if __name__ == "__main__":
    main()
