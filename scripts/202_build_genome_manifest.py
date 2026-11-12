#!/usr/bin/env python3
"""Build ONE verified accession -> genome-path map, and refuse to ship a bad one.

WHY THIS IS ITS OWN SCRIPT. The manifest was previously built inline with
`find ... | xargs cat | sort -u`. Run twice over the same tree it returned
24,353 paths once and 20,947 the next time: xargs splits into several cat
invocations and a failure in any one is silently dropped, so the pipeline's
exit status says nothing. A downstream job then resolved 43% of its genomes to
nothing and still exited 0. A manifest that is silently short is worse than one
that is missing, so this script verifies coverage and exits non-zero rather
than emitting a partial map.

DETERMINISM. Directories are walked with os.scandir and each read is retried,
because the underlying failure is a transient Lustre error, not a missing file.
Every file read is counted and the totals are printed, so a short run is
visible instead of silent.

PRECEDENCE. An accession may appear in several sources. The protctx corpus is
preferred when present because those paths have already been read end to end by
earlier jobs; panel-manifest paths fill the rest. Last-write-wins was a real bug
here: it let an arbitrary duplicate overwrite a known-good path.
"""
import os
import sys
import time

W = "/global/scratch/users/kh36969/fna_ins_discovery"
OUT = W + "/insertions_all"
MAP = OUT + "/genome_map.tsv"


def acc_of(path):
    b = os.path.basename(path)
    for e in (".fna", ".fa", ".fasta"):
        if b.endswith(e):
            return b[:-len(e)]
    return None


def read_retry(path, tries=4):
    for t in range(tries):
        try:
            with open(path) as fh:
                return fh.read()
        except (IOError, OSError) as e:
            if t == tries - 1:
                sys.stderr.write("UNREADABLE %s (%s)\n" % (path, e))
                return None
            time.sleep(0.5 * (t + 1))
    return None


def find_manifests(root):
    """Every panels_manifests/*.manifest, via scandir, with retries."""
    out = []
    try:
        runs = list(os.scandir(root))
    except OSError as e:
        sys.exit("cannot scan %s: %s" % (root, e))
    for d in runs:
        if not d.is_dir():
            continue
        pm = os.path.join(d.path, "panels_manifests")
        if not os.path.isdir(pm):
            continue
        for t in range(4):
            try:
                out.extend(e.path for e in os.scandir(pm)
                           if e.name.endswith(".manifest"))
                break
            except OSError:
                time.sleep(0.5 * (t + 1))
        else:
            sys.exit("could not list %s after retries" % pm)
    return sorted(out)


def main():
    mans = find_manifests(W)
    sys.stderr.write("panel manifest files: %d\n" % len(mans))

    paths, n_read, n_fail = set(), 0, 0
    for i, m in enumerate(mans, 1):
        txt = read_retry(m)
        if txt is None:
            n_fail += 1
            continue
        n_read += 1
        for line in txt.splitlines():
            p = line.strip()
            if p:
                paths.add(p)
        if i % 1000 == 0:
            sys.stderr.write("  %d/%d manifests, %d paths\n" % (i, len(mans), len(paths)))
    sys.stderr.write("manifests read %d, failed %d, distinct paths %d\n"
                     % (n_read, n_fail, len(paths)))
    if n_fail:
        sys.exit("FATAL: %d manifest files unreadable; refusing to emit a "
                 "partial map" % n_fail)

    # protctx corpus takes precedence; panel paths fill the rest
    amap = {}
    for p in sorted(paths):
        a = acc_of(p)
        if a and a not in amap:
            amap[a] = p
    n_panel = len(amap)
    corpus = W + "/protctx/all_fna.txt"
    if os.path.exists(corpus):
        for line in open(corpus):
            p = line.strip()
            a = acc_of(p) if p else None
            if a:
                amap[a] = p                      # preferred source wins
    sys.stderr.write("accessions: %d from panels, %d after corpus merge\n"
                     % (n_panel, len(amap)))

    # VERIFY: every coord_frame the locations table names must resolve
    import sqlite3
    con = sqlite3.connect("file:%s/locations/insertions.sqlite?mode=ro" % W, uri=True)
    want = set(r[0] for r in con.execute(
        "SELECT DISTINCT coord_frame FROM loc WHERE junction_L_abs NOT IN ('.','')"))
    carriers = set(r[0] for r in con.execute(
        "SELECT DISTINCT genome FROM loc WHERE role='carrier_derived'"))
    missing = sorted(want - set(amap))
    miss_car = sorted(carriers - set(amap))
    print("coord_frame genomes wanted : %d   unresolved: %d" % (len(want), len(missing)))
    print("carrier genomes wanted     : %d   unresolved: %d" % (len(carriers), len(miss_car)))
    if missing:
        print("  examples:", missing[:5])

    # VERIFY: the files are actually there
    gone = [a for a in sorted(amap) if not os.path.exists(amap[a])]
    print("accessions whose file is absent on disk: %d" % len(gone))
    if gone:
        print("  examples:", gone[:5])

    with open(MAP, "w") as fh:
        fh.write("accession\tpath\n")
        for a in sorted(amap):
            fh.write("%s\t%s\n" % (a, amap[a]))
    print("wrote %s  (%d accessions)" % (MAP, len(amap)))

    if missing or gone:
        sys.exit("FATAL: manifest does not cover every needed genome")
    print("MANIFEST OK")


if __name__ == "__main__":
    main()
