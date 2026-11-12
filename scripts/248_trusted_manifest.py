#!/usr/bin/env python3
"""248_trusted_manifest.py -- SPEC item 0 (V0 trusted manifest).

Builds, and independently re-checks, the trusted expectations that V2 shard
verification is gated on:

  * the query side: sha256 / byte count / record count of all_queries.fna, the
    sha256 of the sorted query-id list, and the sha256 of all_queries.tsv after
    checking it against the FASTA (V2-3(b) rules);
  * the genome side ("genome trusted manifest", V2-4(c)): for every assembly in
    corpus_v2.tsv its file sha256 and byte count, contig count, base count, the
    sha256 of its contig-id list, and every contig's full id ACC|CONTIG and
    length; plus the 64-way shard partition of V2-4(a).

Subcommands
  build  read the sources, write a NEW run directory, then run `check` on what
         was written.  Exit 0 only if everything passed.
  check  re-read a run directory and verify it against the sources' partition
         and the expected literals, without trusting anything build computed.

Nothing is ever overwritten: build refuses an existing --out directory.  On
failure build writes <out>/FAILED (reason) and never writes v0_result.json;
v0_result.json is written last and only on success.

Parsing rules (must match how C_sweep.sh / V2 turn a genome into shard FASTA):
  * header line: drop the leading '>', split on runs of [ \\t\\n] (awk default
    FS), the first field is the contig id; the derived id is ACC|<contig id>;
  * every other line is sequence; its length is the line without the
    trailing '\\n'.  Sequence lines may contain only ASCII letters.  A '\\r'
    anywhere, a non-letter in a sequence line, sequence before the first
    header, an empty contig id or a zero-length contig is a hard failure
    (we refuse to guess how makeblastdb would treat them).

I/O re-reads (SPEC V0-7; shared scratch has returned EIO / "transport endpoint
shutdown" mid-read):
  * only OSError with errno EIO, ENOTCONN or ESHUTDOWN is retried, and only
    while reading a genome FASTA; every other error fails at once;
  * every read is a complete, fresh read from byte 0 (new hash, new parse);
    nothing from a failed read is kept;
  * at most 3 re-reads per genome (at most 4 reads in all);
  * a genome that saw any such error is accepted only after TWO complete reads
    agree on byte count, sha256 and every (contig id, length).  Each read also
    checks bytes read == file size at open.  The confirming read is one of the
    3 re-reads, so a genome tolerates at most 2 errors, and both must come
    before its first complete read: an error that leaves fewer than 2 reads in
    the budget fails the genome at once (no pointless wait and read);
  * the wait before the re-read that follows error 1 / error 2 is 10 s / 30 s
    (increasing); the confirming read follows a complete read without a wait;
  * errors are counted per node (this process tree runs on one node); at 5 the
    whole build fails and FAILED names the node;
  * every such error is logged the moment it happens (an `IO_ERROR {json}`
    line on stderr, and a shared list the parent reads even after tearing the
    pool down), so a failure record names the genomes, errnos and times;
  * genome_reads.tsv records re-reads and errors for every genome, and
    v0_result.json lists the genomes with re-reads > 0.  Neither is part of the
    genome trusted manifest, so the manifest does not depend on I/O history.
"""
import argparse
import errno
import hashlib
import json
import os
import re
import socket
import subprocess
import sys
import time
from collections import Counter
import multiprocessing

FS_RE = re.compile(rb"[ \t\n]+")
SEQ_OK = re.compile(rb"[A-Za-z]*")
HEX64 = re.compile(r"[0-9a-f]{64}")
QTSV_HEADER = "query_id\tbag_id\tseq_key\tlen_nt"
GENOMES_HEADER = ("shard\tassembly\tfna_path\tfile_bytes\tfile_sha256\t"
                  "n_contigs\tbases\tcontig_ids_sha256")
CONTIGS_HEADER = "shard\tassembly\tfull_id\tlength"
TRUSTED_HEADER = "key\tvalue\tsource"
READS_HEADER = "shard\tassembly\trereads\tio_errors"
# V0-7 I/O re-read policy
IO_RETRY_ERRNOS = {errno.EIO: "EIO", errno.ENOTCONN: "ENOTCONN", errno.ESHUTDOWN: "ESHUTDOWN"}
IO_MAX_REREADS = 3
IO_BACKOFF_SECONDS = (10, 30)           # after error 1, 2 (increasing); see V0-7
NODE_IO_ERROR_LIMIT = 5
INT_KEYS = ("all_queries_fna_bytes", "n_queries", "n_genomes", "n_shards",
            "n_contigs_total", "bases_total")
SHA_KEYS = ("all_queries_fna_sha256", "query_ids_sorted_sha256",
            "all_queries_tsv_sha256", "genome_manifest_sha256",
            "corpus_tsv_sha256")


class Fail(Exception):
    """A check failed; the message names the check."""


class GenomeReadFail(Fail):
    """A genome could not be read; .io holds the re-read record so far."""

    def __init__(self, msg, io):
        super().__init__(msg)
        self.io = io


def sha256_bytes(b):
    return hashlib.sha256(b).hexdigest()


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for blk in iter(lambda: f.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def ids_digest(ids):
    """sha256 of the ids, one per line, each terminated by '\\n'."""
    return sha256_bytes(b"".join(i + b"\n" for i in ids))


def header_id(line, where):
    """Contig/query id from a header line (bytes, '>' included, no '\\n')."""
    toks = [t for t in FS_RE.split(line[1:]) if t]
    if not toks:
        raise Fail(f"empty_id: {where}")
    return toks[0]


def scan_fasta(path, label, opener=open):
    """Read a FASTA file ONCE.  Returns (bytes, sha256, [(id, length)]).

    The byte count read must equal the file size at open time, so a short read
    from the shared filesystem is a failure, not a silently smaller genome."""
    h = hashlib.sha256()
    recs = []
    nread = 0
    cur_id, cur_len = None, 0
    pending = b""
    lineno = 0
    with opener(path, "rb") as f:
        size = os.fstat(f.fileno()).st_size

        def line_done(line):
            nonlocal cur_id, cur_len, lineno
            lineno += 1
            if b"\r" in line:
                raise Fail(f"carriage_return: {label} line {lineno}")
            if line.startswith(b">"):
                if cur_id is not None:
                    if cur_len == 0:
                        raise Fail(f"zero_length_contig: {label} {cur_id.decode(errors='replace')}")
                    recs.append((cur_id, cur_len))
                cur_id, cur_len = header_id(line, f"{label} line {lineno}"), 0
            else:
                if SEQ_OK.fullmatch(line) is None:
                    raise Fail(f"bad_sequence_char: {label} line {lineno}")
                if cur_id is None:
                    if line:
                        raise Fail(f"sequence_before_header: {label} line {lineno}")
                    return
                cur_len += len(line)

        for blk in iter(lambda: f.read(1 << 20), b""):
            nread += len(blk)
            h.update(blk)
            parts = (pending + blk).split(b"\n")
            pending = parts.pop()
            for line in parts:
                line_done(line)
        if pending:
            line_done(pending)
    if nread != size:
        raise Fail(f"short_read: {label} read {nread} of {size} bytes")
    if cur_id is not None:
        if cur_len == 0:
            raise Fail(f"zero_length_contig: {label} {cur_id.decode(errors='replace')}")
        recs.append((cur_id, cur_len))
    if not recs:
        raise Fail(f"no_records: {label}")
    return nread, h.hexdigest(), recs


def read_corpus(path):
    """corpus_v2.tsv -> [(assembly, fna_path)] in file order (header skipped)."""
    rows = []
    with open(path, "rb") as f:
        head = f.readline().rstrip(b"\n").split(b"\t")
        if len(head) < 3 or head[0] != b"assembly" or head[2] != b"fna_path":
            raise Fail("corpus_header")
        for n, line in enumerate(f, start=2):
            if not line.endswith(b"\n"):
                raise Fail(f"corpus_truncated_line: line {n}")
            p = line.rstrip(b"\n").split(b"\t")
            if len(p) < 3 or not p[0] or not p[2]:
                raise Fail(f"corpus_malformed: line {n}")
            rows.append((p[0].decode(), p[2].decode()))
    dups = [a for a, c in Counter(a for a, _ in rows).items() if c > 1]
    if dups:
        raise Fail(f"duplicate_assembly_in_corpus (shards would intersect): {dups[:5]}")
    return rows


def parse_sizes(spec):
    """'35x633,29x634' -> Counter({633: 35, 634: 29})."""
    c = Counter()
    for part in spec.split(","):
        k, _, size = part.partition("x")
        c[int(size)] += int(k)
    return c


def check_partition(rows, n_shards, n_genomes, sizes):
    if len(rows) != n_genomes:
        raise Fail(f"genome_count: {len(rows)} != expected {n_genomes}")
    shards = [set() for _ in range(n_shards)]
    for k, (acc, _) in enumerate(rows):
        shards[k % n_shards].add(acc)
    for i in range(n_shards):
        for j in range(i + 1, n_shards):
            if shards[i] & shards[j]:
                raise Fail(f"shards_intersect: {i} {j}")
    if sum(len(s) for s in shards) != n_genomes or len(set().union(*shards)) != n_genomes:
        raise Fail("shard_union")
    got = Counter(len(s) for s in shards)
    if got != sizes:
        raise Fail(f"shard_sizes: {dict(got)} != expected {dict(sizes)}")


def scan_queries(fna, tsv):
    nbytes, fsha, recs = scan_fasta(fna, "all_queries.fna")
    ids = [i for i, _ in recs]
    dup = [i for i, c in Counter(ids).items() if c > 1]
    if dup:
        raise Fail(f"duplicate_query_id: {[d.decode() for d in dup[:5]]}")
    flen = {i.decode(): n for i, n in recs}
    # all_queries.tsv against the FASTA (V2-3(b))
    seen = set()
    with open(tsv, "rb") as f:
        raw = f.read()
    tsha = sha256_bytes(raw)
    if not raw.endswith(b"\n"):
        raise Fail("queries_tsv_truncated_last_line")
    lines = raw.decode().split("\n")[:-1]
    if not lines or lines[0] != QTSV_HEADER:
        raise Fail("queries_tsv_header")
    for n, line in enumerate(lines[1:], start=2):
        p = line.split("\t")
        if len(p) != 4:
            raise Fail(f"queries_tsv_fields: line {n}")
        qid, bag, key, ln = p
        if qid in seen:
            raise Fail(f"queries_tsv_duplicate_id: {qid}")
        seen.add(qid)
        if qid != "Q_" + key:
            raise Fail(f"queries_tsv_id_not_Q_seqkey: {qid}")
        if not bag:
            raise Fail(f"queries_tsv_empty_bag: {qid}")
        if not re.fullmatch(r"[1-9][0-9]*", ln):
            raise Fail(f"queries_tsv_len_not_int: {qid}")
        if qid not in flen:
            raise Fail(f"queries_tsv_id_not_in_fasta: {qid}")
        if int(ln) != flen[qid]:
            raise Fail(f"queries_tsv_len_mismatch: {qid} {ln} != {flen[qid]}")
    if seen != set(flen):
        raise Fail(f"queries_tsv_id_set != fasta ({len(seen)} vs {len(flen)})")
    return {
        "all_queries_fna_sha256": fsha,
        "all_queries_fna_bytes": nbytes,
        "n_queries": len(ids),
        "query_ids_sorted_sha256": ids_digest(sorted(ids)),
        "all_queries_tsv_sha256": tsha,
    }


def _bump(counter):
    """Add one I/O error to the node counter; return the new node total."""
    with counter.get_lock():
        counter.value += 1
        return counter.value


def _record_io_event(events, label, path, name, attempt, node_total, exc):
    """Log one transient I/O error the moment it happens: a line on stderr (the
    SLURM log) and, when given, an entry in the shared events list, so the
    record survives even if the pool is torn down before the worker returns."""
    ev = {"time": time.strftime("%Y-%m-%dT%H:%M:%S%z"), "host": socket.gethostname(),
          "pid": os.getpid(), "assembly": label, "path": path, "errno": name,
          "attempt": attempt, "node_total": node_total, "message": str(exc)}
    try:
        sys.stderr.write("IO_ERROR " + json.dumps(ev, sort_keys=True) + "\n")
        sys.stderr.flush()
    except OSError:
        pass
    if events is not None:
        try:
            events.append(ev)
        except Exception:   # the record is diagnostic; never let it mask the error itself
            pass


def read_genome(path, label, counter, sleep=time.sleep, events=None):
    """scan_fasta with the V0-7 re-read policy (see the module docstring).

    Returns (nbytes, sha256, recs, io) with io = {"rereads": n, "io_errors": [..]}.
    Raises GenomeReadFail (I/O policy), Fail (content rules) or GenomeReadFail
    wrapping a non-retryable OSError."""
    errs, reads, first = [], 0, None

    def io():
        return {"rereads": max(reads - 1, 0), "io_errors": list(errs)}

    while True:
        if counter.value >= NODE_IO_ERROR_LIMIT:
            raise GenomeReadFail(f"node_io_error_limit_reached: {label} not (fully) read", io())
        reads += 1
        try:
            got = scan_fasta(path, label)
        except FileNotFoundError:
            raise GenomeReadFail(f"genome_missing: {label} {path}", io())
        except OSError as e:
            name = IO_RETRY_ERRNOS.get(e.errno)
            if name is None:
                raise GenomeReadFail(f"genome_unreadable: {label} {path} {e}", io())
            errs.append(name)
            node_total = _bump(counter)
            _record_io_event(events, label, path, name, reads, node_total, e)
            if node_total >= NODE_IO_ERROR_LIMIT:
                raise GenomeReadFail(f"node_io_error_limit: {label} error {name}; node total {node_total}", io())
            # reads left in the budget; after an error two must remain (a
            # complete read and the confirming read)
            if IO_MAX_REREADS + 1 - reads < 2:
                raise GenomeReadFail(f"genome_io_retries_exhausted: {label} {path} errors {','.join(errs)}", io())
            sleep(IO_BACKOFF_SECONDS[len(errs) - 1])
            continue
        if not errs:
            return got + (io(),)
        if first is None:
            # an error happened earlier: a second complete read must agree
            first = got
            continue
        if got != first:
            raise GenomeReadFail(f"genome_reread_mismatch: {label} {path} bytes {first[0]}/{got[0]} "
                                 f"sha256 {first[1]}/{got[1]}", io())
        return got + (io(),)


_NODE_COUNTER = None
_IO_EVENTS = None


def _pool_init(counter, events):
    global _NODE_COUNTER, _IO_EVENTS
    _NODE_COUNTER = counter
    _IO_EVENTS = events


def _scan_genome(job):
    """Pool worker: (k, error or None, (nbytes, sha, recs) or None, io)."""
    k, acc, path = job
    try:
        nbytes, sha, recs, io = read_genome(path, acc, _NODE_COUNTER, events=_IO_EVENTS)
    except GenomeReadFail as e:
        return k, str(e), None, e.io
    except Fail as e:
        return k, str(e), None, None
    ids = [i for i, _ in recs]
    dup = [i for i, c in Counter(ids).items() if c > 1]
    if dup:
        return k, f"duplicate_contig_id: {acc} {[d.decode() for d in dup[:5]]}", None, io
    return k, None, (nbytes, sha, recs), io


def scan_genomes(jobs, threads):
    """Read every genome once (V0-7 re-reads aside), in parallel.

    Returns (res, ios, errors, node_io_errors); res/ios are in job order.  The
    pool is torn down as soon as the node counter reaches the limit."""
    ctx = multiprocessing.get_context("fork")
    counter = ctx.Value("i", 0)
    # a Manager list, not a Queue: workers killed by pool.terminate() cannot
    # corrupt it, and every append is complete when the call returns
    mgr = ctx.Manager()
    events = mgr.list()
    res = [None] * len(jobs)
    ios = [None] * len(jobs)
    errors = []
    with ctx.Pool(threads, initializer=_pool_init, initargs=(counter, events)) as pool:
        for done, (k, err, val, io) in enumerate(pool.imap_unordered(_scan_genome, jobs, chunksize=8), 1):
            ios[k] = io
            if err:
                errors.append(err)
            else:
                res[k] = val
            if counter.value >= NODE_IO_ERROR_LIMIT:
                break
            if done % 5000 == 0:
                print(f"  scanned {done}/{len(jobs)} genomes", flush=True)
    ev = sorted(list(events), key=lambda d: (d["time"], d["assembly"], d["attempt"]))
    mgr.shutdown()
    return res, ios, errors, counter.value, ev


def write_atomic(path, data):
    tmp = path + ".tmp"
    with open(tmp, "wb") as f:
        f.write(data)
        f.flush()
        os.fsync(f.fileno())
    os.rename(tmp, path)


def tool_version(cmd):
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    except OSError as e:
        raise Fail(f"tool_version: {cmd[0]}: {e}")
    if r.returncode != 0 or not r.stdout.strip():
        raise Fail(f"tool_version: {cmd[0]} exit {r.returncode}")
    return r.stdout.strip().splitlines()[0]


def check_blast_version(line, tool, want):
    """G-5 pins BLAST 2.17.0: the measured `<tool> -version` line must name it."""
    m = re.fullmatch(rf"{re.escape(tool)}: ([0-9]+\.[0-9]+\.[0-9]+)\+?", line.strip())
    if not m:
        raise Fail(f"tool_version_unparsable: {tool}: {line!r}")
    if m.group(1) != want:
        raise Fail(f"blast_version_mismatch: {tool} {m.group(1)} != {want} (G-5)")


def git_state(repo, commit, code_files):
    """G-4 git fields, and proof that the code being run is the committed code.

    git_dirty counts untracked files too (an untracked file may be code).  In
    addition every file in code_files must exist in `commit` and the working
    copy must be byte-identical to it; otherwise the run is refused, because
    the recorded commit would not reproduce it."""
    def g(*a, binary=False):
        r = subprocess.run(["git", "-C", repo, *a], capture_output=True, text=not binary)
        if r.returncode != 0:
            err = r.stderr if not binary else r.stderr.decode(errors="replace")
            raise Fail(f"git: {' '.join(a)}: {err.strip()}")
        return r.stdout
    head = g("rev-parse", "HEAD").strip()
    status = g("status", "--porcelain", "--untracked-files=all")
    dirty = bool(status.strip())
    code = {}
    for rel in code_files:
        r = subprocess.run(["git", "-C", repo, "show", f"{commit}:{rel}"], capture_output=True)
        if r.returncode != 0:
            raise Fail(f"code_not_committed: {rel} is not in commit {commit[:12]}")
        committed = sha256_bytes(r.stdout)
        try:
            working = sha256_file(os.path.join(repo, rel))
        except OSError as e:
            raise Fail(f"code_missing: {rel}: {e}")
        if committed != working:
            raise Fail(f"code_differs_from_commit: {rel} (working copy != {commit[:12]})")
        code[rel] = committed
    return head, dirty, status.splitlines(), code


def record_tooltests(d):
    """V0-4: the tool self-tests must have run in this job and passed; record
    every file they wrote (commands, raw output, verdicts) by sha256."""
    need = {"outfmt7_test.verdict": "V0 TEST1 PASS", "v0b_dbtest.log": "V0B PASS",
            "commands.txt": None}
    out = {}
    for name in sorted(os.listdir(d)):
        p = os.path.join(d, name)
        if os.path.isfile(p):
            out[name] = {"bytes": os.path.getsize(p), "sha256": sha256_file(p)}
    for name, token in need.items():
        if name not in out:
            raise Fail(f"tooltest_missing: {name}")
        if token is not None and token not in open(os.path.join(d, name), errors="replace").read():
            raise Fail(f"tooltest_not_passed: {name}")
    with open(os.path.join(d, "commands.txt"), errors="replace") as f:
        out["commands.txt"]["text"] = f.read()
    # the self-test input digests must have been recorded (G-4 / V0-5)
    h = r"[0-9a-f]{64}"
    with open(os.path.join(d, "v0b_dbtest.log"), errors="replace") as f:
        v0b = f.read()
    required = [
        ("commands.txt", out["commands.txt"]["text"], rf"(?m)^# V0 job \S+; test genome \S+ \(sha256 {h}\)$"),
        ("commands.txt", out["commands.txt"]["text"], rf"(?m)^# TEST 1 test inputs: t_q\.fna sha256 {h}, t_db\.fna sha256 {h}$"),
        ("v0b_dbtest.log", v0b, rf"(?m)x\.fna sha256 {h}$"),
        ("v0b_dbtest.log", v0b, rf"(?m)^script sha256: {h}$"),
        ("v0b_dbtest.log", v0b, r"(?m)^git_commit_submit=[0-9a-f]{40} git_commit_runtime=[0-9a-f]{40} git_dirty=(true|false)$"),
    ]
    for name, text, pat in required:
        if not re.search(pat, text):
            raise Fail(f"tooltest_digest_missing: {name}: {pat}")
    return out


def load_expect(path):
    with open(path) as f:
        e = json.load(f)
    need = {"all_queries_fna_sha256", "query_ids_sorted_sha256", "n_queries",
            "all_queries_fna_bytes", "n_genomes", "n_shards", "shard_sizes", "blast_version"}
    miss = need - set(e)
    if miss:
        raise Fail(f"expect_file_missing_keys: {sorted(miss)}")
    extra = set(e) - need
    if extra:
        raise Fail(f"expect_file_unknown_keys: {sorted(extra)}")
    # G-7: counts are integers (JSON int, not float/bool/string), digests are
    # 64 lowercase hex, strings match their grammar; anything else fails
    for k in ("n_queries", "all_queries_fna_bytes", "n_genomes", "n_shards"):
        v = e[k]
        if type(v) is not int or v < 1:
            raise Fail(f"expect_type: {k} must be a positive JSON integer, got {v!r}")
    for k in ("all_queries_fna_sha256", "query_ids_sorted_sha256"):
        if not isinstance(e[k], str) or not HEX64.fullmatch(e[k]):
            raise Fail(f"expect_type: {k} must be 64 lowercase hex, got {e[k]!r}")
    if not isinstance(e["shard_sizes"], str) or not re.fullmatch(r"[1-9][0-9]*x[1-9][0-9]*(,[1-9][0-9]*x[1-9][0-9]*)*", e["shard_sizes"]):
        raise Fail(f"expect_type: shard_sizes must look like '35x633,29x634', got {e['shard_sizes']!r}")
    if not isinstance(e["blast_version"], str) or not re.fullmatch(r"[0-9]+\.[0-9]+\.[0-9]+", e["blast_version"]):
        raise Fail(f"expect_type: blast_version must look like '2.17.0', got {e['blast_version']!r}")
    sizes = parse_sizes(e["shard_sizes"])
    if sum(sizes.values()) != e["n_shards"] or sum(k * v for k, v in sizes.items()) != e["n_genomes"]:
        raise Fail("expect_inconsistent: shard_sizes do not add up to n_shards / n_genomes")
    return e


def compare_literals(q, e):
    for k in ("all_queries_fna_sha256", "query_ids_sorted_sha256",
              "n_queries", "all_queries_fna_bytes"):
        if q[k] != e[k]:
            raise Fail(f"literal_mismatch: {k} got {q[k]} expected {e[k]}")


# ---------------------------------------------------------------- build
def build(a):
    if os.path.lexists(a.out):
        print(f"FAIL out_exists: {a.out} (248 never overwrites a run directory)", file=sys.stderr)
        return 2
    os.makedirs(a.out)
    t0 = time.time()
    host = socket.gethostname()
    try:
        git_submit = os.environ.get("GIT_COMMIT", "")
        if not re.fullmatch(r"[0-9a-f]{40}", git_submit):
            raise Fail("git_commit_submit_unset (pass GIT_COMMIT=$(git rev-parse HEAD) at sbatch)")
        git_runtime, git_dirty, git_status, code_sha = git_state(a.repo, git_submit, a.code_files)
        if os.path.realpath(__file__) != os.path.realpath(os.path.join(a.repo, "scripts/248_trusted_manifest.py")) \
                and "scripts/248_trusted_manifest.py" in a.code_files:
            raise Fail("code_path: the running 248 is not the repo copy")
        tooltests = record_tooltests(a.tooltests)
        tools = {
            "blastn": tool_version([os.path.join(a.blast_bin, "blastn"), "-version"]),
            "makeblastdb": tool_version([os.path.join(a.blast_bin, "makeblastdb"), "-version"]),
            "blastdbcmd": tool_version([os.path.join(a.blast_bin, "blastdbcmd"), "-version"]),
            "python3": sys.version.split()[0],
        }
        e = load_expect(a.expect)
        for t in ("blastn", "makeblastdb", "blastdbcmd"):
            check_blast_version(tools[t], t, e["blast_version"])
        sizes = parse_sizes(e["shard_sizes"])

        # queries
        q = scan_queries(a.queries_fna, a.queries_tsv)
        compare_literals(q, e)
        print(f"queries: {q['n_queries']} records, {q['all_queries_fna_bytes']} bytes; literals match", flush=True)

        # partition
        rows = read_corpus(a.corpus)
        check_partition(rows, e["n_shards"], e["n_genomes"], sizes)
        print(f"partition: {len(rows)} genomes into {e['n_shards']} shards; sizes {dict(sizes)}", flush=True)

        # genomes, read once each (V0-7 re-reads aside), in parallel; results
        # kept in corpus order
        jobs = [(k, acc, path) for k, (acc, path) in enumerate(rows)]
        res, ios, errors, node_io, io_events = scan_genomes(jobs, a.threads)
        ev_short = [f"{d['assembly']}:{d['errno']}@{d['time']}#{d['attempt']}" for d in io_events]
        hit = sorted((rows[k][0], io["rereads"], ",".join(io["io_errors"]))
                     for k, io in enumerate(ios) if io and (io["rereads"] or io["io_errors"]))
        if node_io >= NODE_IO_ERROR_LIMIT:
            raise Fail(f"node_io_error_limit: node {host} had {node_io} I/O errors "
                       f"(limit {NODE_IO_ERROR_LIMIT}); errno counts "
                       f"{dict(Counter(d['errno'] for d in io_events))}; "
                       f"{len(io_events)} logged errors (assembly:errno@time#attempt): {ev_short[:40]}; "
                       f"genomes returned with re-reads: {hit[:20]}")
        if errors:
            errors.sort()
            raise Fail(f"genome_scan: {len(errors)} failures; first: {errors[:5]}; "
                       f"node {host} I/O errors {node_io}; logged errors: {ev_short[:40]}")
        print(f"genomes: {len(rows)} read; node {host} I/O errors {node_io}; "
              f"genomes re-read {len(hit)}: {hit}", flush=True)

        r_lines = [READS_HEADER]
        for k, (acc, _) in enumerate(rows):
            io = ios[k]
            r_lines.append(f"{k % e['n_shards']}\t{acc}\t{io['rereads']}\t"
                           f"{','.join(io['io_errors']) or '-'}")
        g_lines = [GENOMES_HEADER]
        c_lines = [CONTIGS_HEADER]
        shard_ids = [set() for _ in range(e["n_shards"])]
        n_contigs_total = bases_total = 0
        for k, (acc, path) in enumerate(rows):
            nbytes, sha, recs = res[k]
            s = k % e["n_shards"]
            bases = sum(n for _, n in recs)
            n_contigs_total += len(recs)
            bases_total += bases
            g_lines.append(f"{s}\t{acc}\t{path}\t{nbytes}\t{sha}\t{len(recs)}\t{bases}\t"
                           f"{ids_digest([i for i, _ in recs])}")
            for cid, n in recs:
                full = f"{acc}|{cid.decode()}"
                if full in shard_ids[s]:
                    raise Fail(f"duplicate_full_id_in_shard: shard {s} {full}")
                shard_ids[s].add(full)
                c_lines.append(f"{s}\t{acc}\t{full}\t{n}")
        g_bytes = ("\n".join(g_lines) + "\n").encode()
        c_bytes = ("\n".join(c_lines) + "\n").encode()
        write_atomic(os.path.join(a.out, "genomes.tsv"), g_bytes)
        write_atomic(os.path.join(a.out, "contigs.tsv"), c_bytes)
        index = f"{sha256_bytes(g_bytes)}  genomes.tsv\n{sha256_bytes(c_bytes)}  contigs.tsv\n".encode()
        write_atomic(os.path.join(a.out, "genome_manifest.index"), index)
        write_atomic(os.path.join(a.out, "genome_reads.tsv"), ("\n".join(r_lines) + "\n").encode())

        src = f"V0 run {os.environ.get('SLURM_JOB_ID', 'local')} (248 build, git {git_submit[:12]})"
        lit = "实测：V0 job 27769455；本次 V0 重新计算后逐值相等"
        trusted = [
            ("all_queries_fna_sha256", q["all_queries_fna_sha256"], lit),
            ("query_ids_sorted_sha256", q["query_ids_sorted_sha256"], lit),
            ("n_queries", q["n_queries"], lit),
            ("all_queries_fna_bytes", q["all_queries_fna_bytes"], "实测：job 27768584 暂存时；本次 V0 重新计算后相等"),
            ("all_queries_tsv_sha256", q["all_queries_tsv_sha256"], src),
            ("genome_manifest_sha256", sha256_bytes(index), src + "；= sha256(genome_manifest.index)"),
            ("corpus_tsv_sha256", sha256_file(a.corpus), src),
            ("n_genomes", len(rows), src),
            ("n_shards", e["n_shards"], src),
            ("n_contigs_total", n_contigs_total, src),
            ("bases_total", bases_total, src),
        ]
        t_bytes = (TRUSTED_HEADER + "\n" + "".join(f"{k}\t{v}\t{s}\n" for k, v, s in trusted)).encode()
        write_atomic(os.path.join(a.out, "isref_all_trusted.tsv"), t_bytes)

        # independent re-check of what is now on disk
        check_dir(a.out, a.corpus, e)
        print("check: run directory re-read and verified", flush=True)

        files = {}
        for name in RESULT_OUTPUTS:
            p = os.path.join(a.out, name)
            files[name] = {"bytes": os.path.getsize(p), "sha256": sha256_file(p)}
        result = {
            "verified": True,
            "git_commit_submit": git_submit,
            "git_commit_runtime": git_runtime,
            "git_dirty": git_dirty,
            "git_status_porcelain": git_status,
            "code_files_sha256_at_commit": code_sha,
            "tooltests": tooltests,
            "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
            "slurm_array_task_id": os.environ.get("SLURM_ARRAY_TASK_ID"),
            "host": host,
            "tools": tools,
            "io_rereads": {
                "policy": {"retry_errnos": sorted(IO_RETRY_ERRNOS.values()),
                           "max_rereads_per_genome": IO_MAX_REREADS,
                           "backoff_seconds": list(IO_BACKOFF_SECONDS),
                           "node_io_error_limit": NODE_IO_ERROR_LIMIT,
                           "confirming_read_after_error": True},
                "node": host,
                "node_io_errors": node_io,
                "events": io_events,
                "genomes_reread": [{"assembly": acc, "rereads": n, "io_errors": errs}
                                   for acc, n, errs in hit],
            },
            "argv": sys.argv,
            "expect": e,
            "inputs": {
                "all_queries.fna": {"path": a.queries_fna, "sha256": q["all_queries_fna_sha256"],
                                    "bytes": q["all_queries_fna_bytes"]},
                "all_queries.tsv": {"path": a.queries_tsv, "sha256": q["all_queries_tsv_sha256"]},
                "corpus_v2.tsv": {"path": a.corpus, "sha256": sha256_file(a.corpus)},
                "expect": {"path": a.expect, "sha256": sha256_file(a.expect)},
                "script_248_sha256": sha256_file(os.path.abspath(__file__)),
            },
            "counts": {"queries": q["n_queries"], "genomes": len(rows),
                       "contigs": n_contigs_total, "bases": bases_total},
            "outputs": files,
            "seconds": round(time.time() - t0, 1),
        }
        write_atomic(os.path.join(a.out, "v0_result.json"),
                     (json.dumps(result, indent=1, sort_keys=True) + "\n").encode())
        print(f"V0 BUILD VERIFIED in {result['seconds']}s -> {a.out}")
        return 0
    except (Fail, OSError) as ex:
        msg = str(ex) if isinstance(ex, Fail) else f"os_error: {ex!r}"
        write_atomic(os.path.join(a.out, "FAILED"),
                     f"{msg}\nhost={host} slurm_job_id={os.environ.get('SLURM_JOB_ID')}\n".encode())
        print(f"FAIL {msg} (host {host})", file=sys.stderr)
        return 1


# ---------------------------------------------------------------- check
def _int(v, what):
    if not re.fullmatch(r"(0|[1-9][0-9]*)", v):
        raise Fail(f"type: {what} not a decimal integer: {v!r}")
    return int(v)


def _sha(v, what):
    if not HEX64.fullmatch(v):
        raise Fail(f"type: {what} not 64 lowercase hex: {v!r}")
    return v


def _read_lines(path):
    with open(path, "rb") as f:
        raw = f.read()
    if not raw.endswith(b"\n"):
        raise Fail(f"truncated_last_line: {os.path.basename(path)}")
    return raw, raw.decode().split("\n")[:-1]


RESULT_OUTPUTS = ("genomes.tsv", "contigs.tsv", "genome_manifest.index", "isref_all_trusted.tsv",
                  "genome_reads.tsv")


def check_dir(out, corpus, e, require_result=False):
    """Verify a run directory from scratch.  Raises Fail.

    require_result=True (the `check` subcommand, i.e. a finished run): also
    require v0_result.json and verify that every output it records still has
    the recorded bytes and sha256, and that its io_rereads summary equals what
    genome_reads.tsv says.  build calls this before v0_result.json exists."""
    rows = read_corpus(corpus)
    n_shards = e["n_shards"]
    check_partition(rows, n_shards, e["n_genomes"], parse_sizes(e["shard_sizes"]))

    idx_raw, idx = _read_lines(os.path.join(out, "genome_manifest.index"))
    want = {}
    for line in idx:
        m = re.fullmatch(r"([0-9a-f]{64})  (genomes\.tsv|contigs\.tsv)", line)
        if not m or m.group(2) in want:
            raise Fail(f"index_malformed: {line!r}")
        want[m.group(2)] = m.group(1)
    if set(want) != {"genomes.tsv", "contigs.tsv"}:
        raise Fail("index_incomplete")
    g_raw, g = _read_lines(os.path.join(out, "genomes.tsv"))
    c_raw, c = _read_lines(os.path.join(out, "contigs.tsv"))
    if sha256_bytes(g_raw) != want["genomes.tsv"]:
        raise Fail("index_sha_mismatch: genomes.tsv")
    if sha256_bytes(c_raw) != want["contigs.tsv"]:
        raise Fail("index_sha_mismatch: contigs.tsv")

    if not g or g[0] != GENOMES_HEADER:
        raise Fail("genomes_header")
    if len(g) - 1 != len(rows):
        raise Fail(f"genomes_rows: {len(g) - 1} != corpus {len(rows)}")
    genomes = []
    for k, line in enumerate(g[1:]):
        p = line.split("\t")
        if len(p) != 8:
            raise Fail(f"genomes_fields: row {k}")
        s = _int(p[0], "shard")
        if (p[1], p[2]) != rows[k] or s != k % n_shards:
            raise Fail(f"genomes_row_vs_corpus: row {k} {p[1]}")
        _int(p[3], "file_bytes")
        _sha(p[4], "file_sha256")
        nc = _int(p[5], "n_contigs")
        nb = _int(p[6], "bases")
        _sha(p[7], "contig_ids_sha256")
        if nc < 1 or nb < 1:
            raise Fail(f"genomes_empty: {p[1]}")
        genomes.append((s, p[1], nc, nb, p[7]))

    if not c or c[0] != CONTIGS_HEADER:
        raise Fail("contigs_header")
    pos = 1
    shard_ids = [set() for _ in range(n_shards)]
    n_contigs_total = bases_total = 0
    for s, acc, nc, nb, idsha in genomes:
        block = c[pos:pos + nc]
        pos += nc
        if len(block) != nc:
            raise Fail(f"contig_count: {acc} has {len(block)} rows, genomes.tsv says {nc}")
        ids, total = [], 0
        for line in block:
            p = line.split("\t")
            if len(p) != 4:
                raise Fail(f"contigs_fields: {acc}")
            if _int(p[0], "shard") != s or p[1] != acc:
                raise Fail(f"contig_row_out_of_order: expected {acc}, got {p[1]}")
            if not p[2].startswith(acc + "|") or len(p[2]) == len(acc) + 1:
                raise Fail(f"contig_full_id: {p[2]}")
            if p[2] in shard_ids[s]:
                raise Fail(f"duplicate_full_id_in_shard: shard {s} {p[2]}")
            shard_ids[s].add(p[2])
            n = _int(p[3], "length")
            if n < 1:
                raise Fail(f"zero_length_contig: {p[2]}")
            total += n
            ids.append(p[2][len(acc) + 1:].encode())
        if total != nb:
            raise Fail(f"contig_length_sum: {acc} sum {total} != bases {nb}")
        if ids_digest(ids) != idsha:
            raise Fail(f"contig_ids_sha: {acc}")
        n_contigs_total += nc
        bases_total += nb
    if pos != len(c):
        raise Fail(f"contigs_extra_rows: {len(c) - pos}")

    # genome_reads.tsv (V0-7): one row per genome in corpus order; re-reads are
    # 0 without an error, and errors + 1 (the confirming read) with errors
    _, r = _read_lines(os.path.join(out, "genome_reads.tsv"))
    if not r or r[0] != READS_HEADER or len(r) - 1 != len(rows):
        raise Fail("genome_reads_header_or_rows")
    node_io, reread = 0, []
    for k, line in enumerate(r[1:]):
        p = line.split("\t")
        if len(p) != 4 or p[1] != rows[k][0] or _int(p[0], "shard") != k % n_shards:
            raise Fail(f"genome_reads_row: {k}")
        n = _int(p[2], "rereads")
        errs = [] if p[3] == "-" else p[3].split(",")
        if any(x not in IO_RETRY_ERRNOS.values() for x in errs):
            raise Fail(f"genome_reads_errno: {p[1]} {p[3]!r}")
        if n != (len(errs) + 1 if errs else 0) or n > IO_MAX_REREADS:
            raise Fail(f"genome_reads_count: {p[1]} rereads {n} errors {p[3]}")
        node_io += len(errs)
        if n:
            reread.append({"assembly": p[1], "rereads": n, "io_errors": p[3]})
    if node_io >= NODE_IO_ERROR_LIMIT:
        raise Fail(f"genome_reads_node_limit: {node_io}")

    t_raw, t = _read_lines(os.path.join(out, "isref_all_trusted.tsv"))
    if not t or t[0] != TRUSTED_HEADER:
        raise Fail("trusted_header")
    kv = {}
    for line in t[1:]:
        p = line.split("\t")
        if len(p) != 3 or not p[2] or p[0] in kv:
            raise Fail(f"trusted_row: {line!r}")
        kv[p[0]] = p[1]
    if set(kv) != set(INT_KEYS) | set(SHA_KEYS):
        raise Fail(f"trusted_keys: {sorted(set(kv) ^ (set(INT_KEYS) | set(SHA_KEYS)))}")
    v = {k: _int(kv[k], k) for k in INT_KEYS}
    v.update({k: _sha(kv[k], k) for k in SHA_KEYS})
    for k in ("all_queries_fna_sha256", "query_ids_sorted_sha256", "n_queries", "all_queries_fna_bytes"):
        if v[k] != e[k]:
            raise Fail(f"literal_mismatch: {k} got {v[k]} expected {e[k]}")
    if v["genome_manifest_sha256"] != sha256_bytes(idx_raw):
        raise Fail("trusted_genome_manifest_sha256 != sha256(genome_manifest.index)")
    if v["corpus_tsv_sha256"] != sha256_file(corpus):
        raise Fail("trusted_corpus_sha256")
    if (v["n_genomes"], v["n_shards"], v["n_contigs_total"], v["bases_total"]) != \
            (len(rows), n_shards, n_contigs_total, bases_total):
        raise Fail("trusted_totals")

    if not require_result:
        return
    rp = os.path.join(out, "v0_result.json")
    if not os.path.exists(rp):
        raise Fail("result_missing: v0_result.json (run incomplete or failed)")
    try:
        with open(rp, "rb") as f:
            res = json.loads(f.read())
    except ValueError as ex:
        raise Fail(f"result_unparsable: {ex}")
    if not isinstance(res, dict) or res.get("verified") is not True:
        raise Fail("result_not_verified")
    outs = res.get("outputs")
    if not isinstance(outs, dict) or set(outs) != set(RESULT_OUTPUTS):
        raise Fail(f"result_outputs_keys: {sorted(outs) if isinstance(outs, dict) else outs!r}")
    for name in RESULT_OUTPUTS:
        p = os.path.join(out, name)
        rec = outs[name]
        if not isinstance(rec, dict) or rec.get("bytes") != os.path.getsize(p) \
                or rec.get("sha256") != sha256_file(p):
            raise Fail(f"result_output_digest: {name} differs from v0_result.json")
    io = res.get("io_rereads")
    if not isinstance(io, dict) or io.get("node_io_errors") != node_io \
            or io.get("genomes_reread") != sorted(reread, key=lambda d: d["assembly"]):
        raise Fail("result_io_rereads != genome_reads.tsv")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build")
    b.add_argument("--queries-fna", required=True)
    b.add_argument("--queries-tsv", required=True)
    b.add_argument("--corpus", required=True)
    b.add_argument("--expect", required=True, help="JSON with the SPEC V2-2 literals and partition sizes")
    b.add_argument("--blast-bin", required=True)
    b.add_argument("--repo", required=True)
    b.add_argument("--code-files", nargs="+", required=True,
                   help="repo-relative files that must be byte-identical to GIT_COMMIT")
    b.add_argument("--tooltests", required=True, help="directory with the V0-4 tool self-test outputs")
    b.add_argument("--threads", type=int, default=8)
    b.add_argument("--out", required=True)
    c = sub.add_parser("check")
    c.add_argument("--out", required=True)
    c.add_argument("--corpus", required=True)
    c.add_argument("--expect", required=True)
    a = ap.parse_args()
    if a.cmd == "build":
        return build(a)
    try:
        check_dir(a.out, a.corpus, load_expect(a.expect), require_result=True)
    except Fail as ex:
        print(f"FAIL {ex}", file=sys.stderr)
        return 1
    print(f"CHECK PASS {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
