#!/usr/bin/env python3
"""249_io_probe.py -- read-only I/O probe of the genome FASTA files on Lustre.

Diagnostic only (not a SPEC item; no output is used by the pipeline).  Run by
slurm/isref_all/IO_probe.sh on one compute node.  Reads the corpus genome files
and writes ONLY into --out (a new directory).

  B  reads: disjoint random samples (fixed seed) read in full at concurrency
     1, 8, 32, 64; every file's open/read time, bytes and any error (errno,
     wall-clock time) is logged as it happens.  Each stage has a wall budget;
     after it no new files are started.
  C  layout: `lfs getstripe -y` of every sampled file -> the OST indices it
     uses (PFL: first 64 KiB on the NVMe pool, the rest on HDD OSTs)
  D  write: 512 MiB written with fsync into --out, then read back
  A  metadata (last, with its own wall budget): os.stat every corpus fna_path
     (follows symlinks), 32 threads
Summary: throughput per stage, latency percentiles, errors by errno, per-OST
read speed and errors, per-directory stat latency and errors.
"""
import argparse
import errno
import json
import os
import random
import re
import statistics
import subprocess
import sys
import threading
import time
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor

BUF = 8 << 20
LOCK = threading.Lock()


def now():
    return time.strftime("%H:%M:%S")


def log(msg):
    with LOCK:
        print(f"[{now()}] {msg}", flush=True)


def ename(e):
    return errno.errorcode.get(e.errno, str(e.errno)) if isinstance(e, OSError) else type(e).__name__


def pct(xs, p):
    if not xs:
        return None
    xs = sorted(xs)
    return round(xs[min(len(xs) - 1, int(p / 100 * len(xs)))], 3)


STAT_DEADLINE = [float("inf")]


def stat_one(job):
    k, acc, path = job
    if time.time() > STAT_DEADLINE[0]:
        return k, acc, path, None, None, "SKIPPED_BUDGET", None
    t = time.time()
    try:
        st = os.stat(path)
        return k, acc, path, round(time.time() - t, 4), st.st_size, None, None
    except OSError as e:
        log(f"STAT_ERROR {acc} {path} {ename(e)} {e}")
        return k, acc, path, round(time.time() - t, 4), None, ename(e), now()


def read_one(job):
    stage, acc, path = job
    t0 = time.time()
    n = 0
    t_open = None
    try:
        with open(path, "rb", buffering=0) as f:
            t_open = time.time() - t0
            while True:
                b = f.read(BUF)
                if not b:
                    break
                n += len(b)
        dt = time.time() - t0
        return dict(stage=stage, acc=acc, path=path, bytes=n, open_s=round(t_open, 4),
                    total_s=round(dt, 4), err=None, err_time=None)
    except OSError as e:
        dt = time.time() - t0
        log(f"READ_ERROR stage={stage} {acc} {path} {ename(e)} after {n} bytes {dt:.1f}s: {e}")
        return dict(stage=stage, acc=acc, path=path, bytes=n,
                    open_s=None if t_open is None else round(t_open, 4),
                    total_s=round(dt, 4), err=ename(e), err_time=now())


def run_stage(stage, conc, files, budget_s):
    """Read files at concurrency conc; stop starting new files after budget_s."""
    log(f"stage {stage}: concurrency {conc}, {len(files)} files, budget {budget_s}s")
    out = []
    deadline = time.time() + budget_s
    t0 = time.time()
    it = iter(files)
    pending = set()
    with ThreadPoolExecutor(conc) as ex:
        for _ in range(conc):
            nxt = next(it, None)
            if nxt is None:
                break
            pending.add(ex.submit(read_one, (stage, *nxt)))
        while pending:
            done = {f for f in pending if f.done()}
            if not done:
                time.sleep(0.05)
                continue
            for f in done:
                pending.remove(f)
                out.append(f.result())
                if time.time() < deadline:
                    nxt = next(it, None)
                    if nxt is not None:
                        pending.add(ex.submit(read_one, (stage, *nxt)))
            if len(out) % 200 == 0 and done:
                log(f"  stage {stage}: {len(out)} done")
    wall = time.time() - t0
    ok = [r for r in out if r["err"] is None]
    tot = sum(r["bytes"] for r in out)
    s = dict(stage=stage, concurrency=conc, files_planned=len(files), files_done=len(out),
             files_ok=len(ok), errors=dict(Counter(r["err"] for r in out if r["err"])),
             wall_s=round(wall, 1), bytes=tot, MB_per_s=round(tot / 1e6 / wall, 1) if wall else None,
             per_file_MB_per_s_p50=pct([r["bytes"] / 1e6 / r["total_s"] for r in ok if r["total_s"] > 0], 50),
             open_s_p50=pct([r["open_s"] for r in ok], 50), open_s_p99=pct([r["open_s"] for r in ok], 99),
             total_s_p50=pct([r["total_s"] for r in ok], 50), total_s_p90=pct([r["total_s"] for r in ok], 90),
             total_s_p99=pct([r["total_s"] for r in ok], 99), total_s_max=pct([r["total_s"] for r in ok], 100),
             budget_hit=len(out) < len(files))
    log(f"stage {stage} summary: {json.dumps(s)}")
    return out, s


def osts(path):
    try:
        r = subprocess.run(["lfs", "getstripe", "-y", os.path.realpath(path)],
                           capture_output=True, text=True, timeout=60)
        return sorted({int(x) for x in re.findall(r"l_ost_idx:\s*(\d+)", r.stdout)})
    except Exception as e:  # diagnostic: record and continue
        return [f"ERR:{type(e).__name__}"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--corpus", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--stage-budget", type=int, default=420)
    ap.add_argument("--stat-budget", type=int, default=360)
    a = ap.parse_args()
    os.makedirs(a.out)          # must be new
    rows = []
    with open(a.corpus) as f:
        f.readline()
        for k, line in enumerate(f):
            p = line.rstrip("\n").split("\t")
            rows.append((k, p[0], p[2]))
    log(f"corpus: {len(rows)} genomes")
    summary = {"host": os.uname().nodename, "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
               "git_commit_submit": os.environ.get("GIT_COMMIT"), "n_genomes": len(rows)}

    # B: reads at increasing concurrency, disjoint samples
    ok_rows = [(acc, path) for k, acc, path in rows]
    random.Random(a.seed).shuffle(ok_rows)
    plan = [("c1", 1, 40), ("c8", 8, 320), ("c32", 32, 960), ("c64", 64, 1920)]
    reads, pos, stages = [], 0, []
    for stage, conc, n in plan:
        out, s = run_stage(stage, conc, ok_rows[pos:pos + n], a.stage_budget)
        pos += n
        reads += out
        stages.append(s)
        summary["B_reads"] = stages
        with open(os.path.join(a.out, "summary.partial.json"), "w") as f:
            json.dump(summary, f, indent=1)

    # C: layout of every sampled file
    log(f"C: lfs getstripe for {len(reads)} files")
    with ThreadPoolExecutor(16) as ex:
        lays = list(ex.map(lambda r: osts(r["path"]), reads))
    with open(os.path.join(a.out, "reads.tsv"), "w") as f:
        f.write("stage\tassembly\tpath\tbytes\topen_s\ttotal_s\terr\terr_time\tosts\n")
        for r, o in zip(reads, lays):
            r["osts"] = o
            f.write("\t".join("" if r[c] is None else str(r[c]) for c in
                              ("stage", "acc", "path", "bytes", "open_s", "total_s", "err", "err_time"))
                    + "\t" + ",".join(map(str, o)) + "\n")
    per = defaultdict(lambda: {"n": 0, "err": 0, "mbps": []})
    for r in reads:
        hdd = [o for o in r["osts"] if isinstance(o, int) and o < 100]
        key = f"OST{hdd[0]}" if len(hdd) == 1 else ("nvme_only" if not hdd else "multi:" + ",".join(map(str, hdd)))
        per[key]["n"] += 1
        if r["err"]:
            per[key]["err"] += 1
        elif r["total_s"] > 0 and r["bytes"] > 1e6:
            per[key]["mbps"].append(r["bytes"] / 1e6 / r["total_s"])
    summary["C_per_hdd_ost"] = {k: {"n": v["n"], "errors": v["err"],
                                    "file_MB_per_s_p50": pct(v["mbps"], 50),
                                    "file_MB_per_s_p10": pct(v["mbps"], 10)}
                                for k, v in sorted(per.items())}

    # D: write test
    p = os.path.join(a.out, "write_test.bin")
    blk = os.urandom(8 << 20)
    try:
        t0 = time.time()
        with open(p, "wb") as f:
            for _ in range(64):
                f.write(blk)
            f.flush()
            os.fsync(f.fileno())
        tw = time.time() - t0
        t0 = time.time()
        n = 0
        with open(p, "rb") as f:
            while True:
                b = f.read(BUF)
                if not b:
                    break
                n += len(b)
        tr = time.time() - t0
        summary["D_write"] = {"bytes": 64 * len(blk), "write_fsync_s": round(tw, 2),
                              "write_MB_per_s": round(64 * len(blk) / 1e6 / tw, 1),
                              "readback_s": round(tr, 2), "readback_ok": n == 64 * len(blk)}
    except OSError as e:
        summary["D_write"] = {"error": ename(e), "msg": str(e), "time": now()}
    finally:
        try:
            os.remove(p)
        except OSError:
            pass
    log(f"D write summary: {json.dumps(summary['D_write'])}")

    # A: metadata
    t0 = time.time()
    STAT_DEADLINE[0] = t0 + a.stat_budget
    with ThreadPoolExecutor(32) as ex:
        meta = list(ex.map(stat_one, rows, chunksize=64))
    a_wall = time.time() - t0
    with open(os.path.join(a.out, "stat.tsv"), "w") as f:
        f.write("k\tassembly\tpath\tstat_s\tsize\terr\terr_time\n")
        for m in meta:
            f.write("\t".join("" if x is None else str(x) for x in m) + "\n")
    bydir = defaultdict(list)
    for m in meta:
        bydir[os.path.dirname(m[2])].append(m)
    summary["A_stat"] = {
        "wall_s": round(a_wall, 1), "errors": dict(Counter(m[5] for m in meta if m[5])),
        "stat_s_p50": pct([m[3] for m in meta if m[3] is not None], 50),
        "stat_s_p99": pct([m[3] for m in meta if m[3] is not None], 99),
        "stat_s_max": pct([m[3] for m in meta if m[3] is not None], 100),
        "total_bytes": sum(m[4] for m in meta if m[4]),
        "by_dir": {d: {"n": len(v), "errors": sum(1 for m in v if m[5] and m[5] != "SKIPPED_BUDGET"),
                       "stat_s_p50": pct([m[3] for m in v if m[3] is not None], 50),
                       "stat_s_max": pct([m[3] for m in v if m[3] is not None], 100)}
                   for d, v in sorted(bydir.items(), key=lambda kv: -len(kv[1]))}}
    log(f"A stat summary: {json.dumps({k: v for k, v in summary['A_stat'].items() if k != 'by_dir'})}")

    with open(os.path.join(a.out, "summary.json"), "w") as f:
        json.dump(summary, f, indent=1)
    log("PROBE DONE")


if __name__ == "__main__":
    sys.exit(main())
