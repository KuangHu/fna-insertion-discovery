#!/usr/bin/env python3
"""ARM B -- cross-genome empty-vs-filled allele discovery. Family-agnostic.

Builds a sequence-based pangenome graph over ONE panel of homologous
assemblies, then reads the bubbles. A bubble whose alleles differ in length by
>= --min-insert IS the empty/filled structure:

    genome A:  LEFT ------------- RIGHT     short allele (empty)
    genome B:  LEFT -- INSERT --- RIGHT     long allele  (filled)

Nothing here knows what an IS is. No transposase model, no IS library, no
annotation. Output is a STRUCTURAL ALLELE, not an insertion: with assemblies
alone gain cannot be told from loss without an outgroup.

Engine: minigraph (-cxggs graph, -cxasm --call genotyping) + gfatools bubble.

SPEC item M1-B (docs/SPEC.md) is the contract this script is held to:

  * B-1  One run per output directory (exclusive flock on <outdir>/.armB.lock);
         the old completion marker and outputs are removed before any check
         that can fail. Every output is written into a per-run
         <outdir>/.armB_tmp.<random>/ and moved to its
         final name only after the whole run succeeded; armB_complete.json is
         written LAST. Stale final outputs are removed at start, so a crash
         never leaves a non-empty final file. Any tool failure, bad input or
         unparseable line exits non-zero with a FATAL message. --reuse reuses
         the graph only if it was built from the same inputs, and a call file
         only if it was made against the same graph and the same FASTA, both
         identified by content (sha256); and each only if the minigraph that
         made it is the one running now (measured version and executable
         sha256) with the same arguments. The marker says, for the graph and
         every call file, whether it was reused and which run created it
         (generator, time, commit, host). Every call file must hold exactly one
         record per bubble; an absent record is a truncated file, not a
         missing call. GIT_COMMIT (40 hex) is required (G-4).
  * B-2  `gfatools bubble` lines must have the 14 documented columns and
         `minigraph --call` records the 6 documented fields; anything else is
         a parse failure (exit 2), never a guess. A missing call (`*`, `.`) is
         kept apart from a zero-length allele.
  * B-3  Every genome at a bubble is empty, filled, intermediate or missing;
         all four counts are written. Thresholds live in classify() only.
  * B-4  min_insert <= insert_size <= max_insert -> armB_events.tsv; larger ->
         armB_large_events.tsv with its sequence; smaller -> counted in the
         summary. Nothing is dropped silently.
  * B-5  Allele sequences come from the CALLED genomes' own graph walks (GFA
         segments), not from the bubble's shortest/longest path. The exported
         insert is diff(empty reference walk, filled reference walk); a
         length mismatch with insert_size beyond --len-tol, and bubbles with
         more than two length classes, are flagged.
  * B-6  Junctions are 0-based half-open contig coordinates of
         empty_ref_genome, given ONLY when that genome's own sequence was
         located exactly and uniquely: source-node tail + its walk + sink-node
         head must occur exactly once in the whole contig (both strands), on
         the call's strand, overlapping the call interval. Otherwise junction_L/R = -1,
         frame_status says why, and the raw minigraph call interval is given
         in separate columns (it is approximate: measured slack up to ~570 bp).
  * B-7  TSD is annotation only; no code path filters on it.
  * B-9  0 bubbles exits 3 unless --allow-zero-bubbles.

Outputs (in --outdir): graph.gfa (+ .meta.json), calls/<sample>.bed (+ .meta.json),
bubbles.tsv, armB_events.tsv, armB_pairs.tsv, armB_inserts.fna,
armB_empty_sites.fna, armB_large_events.tsv, armB_large_inserts.fna,
armB_failed.json (instead of the marker, on any failure: reason, counts such
as unparseable lines, G-4 context), armB_complete.json (written last; G-4: git commits, host, SLURM ids, tool
versions measured at run time, parameters, input sha256, counts, output sha256).
"""
import argparse
import fcntl
import hashlib
import json
import os
import re
import shutil
import statistics
import subprocess
import sys
import tempfile
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib.fasta import Fasta, revcomp                            # noqa: E402
from lib.tsd import diff_alleles, find_tsd                      # noqa: E402
from lib.util import log, run, seq_md5                          # noqa: E402

MISSING = {".", "*"}          # an EMPTY field is a truncated record, not a missing call
WALK_TOKEN = re.compile(r"([<>])([^<>]+)")
TMPDIR_PREFIX = ".armB_tmp."
LOCKFILE = ".armB.lock"
SEGLIST = re.compile(r"[^,\s]+(,[^,\s]+)*")
SEQ = re.compile(r"[A-Za-z]+")
MARKER = "armB_complete.json"
FAILED = "armB_failed.json"
OUTPUTS = ("bubbles.tsv", "armB_events.tsv", "armB_pairs.tsv", "armB_inserts.fna",
           "armB_empty_sites.fna", "armB_large_events.tsv", "armB_large_inserts.fna")
EVENT_COLS = [
    "event_id", "cluster_id", "chrom", "bubble_start", "bubble_end",
    "is_inversion", "len_empty", "len_filled", "insert_size",
    "n_empty", "n_filled", "n_intermediate", "n_missing", "n_length_classes",
    "allele_state", "empty_ref_genome", "empty_ref_contig", "empty_ref_strand",
    "empty_ref_call_start", "empty_ref_call_end", "frame_status",
    "junction_L", "junction_R", "junction_span", "architecture",
    "filled_ref_genome", "insert_len_direct", "insert_md5",
    "tsd_len", "tsd_seq", "tsd_side", "tsd_conf",
    "empty_genomes", "filled_genomes", "qc_tier", "flags"]
LARGE_COLS = [
    "large_event_id", "cluster_id", "chrom", "bubble_start", "bubble_end",
    "is_inversion", "len_empty", "len_filled", "insert_size",
    "n_empty", "n_filled", "n_intermediate", "n_missing", "n_length_classes",
    "empty_ref_genome", "filled_ref_genome", "insert_len_direct", "insert_md5",
    "empty_genomes", "filled_genomes", "flags"]


class ArmBFail(Exception):
    """A FATAL condition: message for the log, exit code for the caller."""

    def __init__(self, msg, code=2, details=None):
        super().__init__(msg)
        self.code = code
        self.details = details or {}


def sample_name(path):
    b = os.path.basename(path)
    for suf in (".fna.gz", ".fa.gz", ".fasta.gz", ".fna", ".fa", ".fasta"):
        if b.endswith(suf):
            return b[:-len(suf)]
    return b


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for blk in iter(lambda: fh.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


_IDENTITY = {}


def file_identity(path):
    """Content identity for --reuse (B-1): real path, size and sha256. Cached
    per (path, size, mtime) within one run so each input is hashed once."""
    st = os.stat(path)
    key = (os.path.realpath(path), st.st_size, st.st_mtime_ns)
    if key not in _IDENTITY:
        _IDENTITY[key] = [key[0], st.st_size, sha256_file(path)]
    return _IDENTITY[key]


def read_json(path):
    try:
        with open(path) as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return None


def write_json(path, obj):
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh, indent=1, sort_keys=True)
        fh.write("\n")
    os.replace(tmp, path)


def tool_versions(minigraph, gfatools):
    """G-4: measured at run time, never hard-coded; a failure is fatal."""
    v = {"python3": sys.version.split()[0]}
    for name, cmd in (("minigraph", [minigraph, "--version"]), ("gfatools", [gfatools, "version"])):
        out = run(cmd).stdout.strip().splitlines()
        if not out:
            raise ArmBFail("%s printed no version" % name)
        v[name] = out[-1]
    return v


def run_context():
    return {"host": os.uname().nodename,
            "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
            "slurm_array_task_id": os.environ.get("SLURM_ARRAY_TASK_ID")}


def git_provenance():
    repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    def g(*a):
        try:
            return subprocess.run(["git", "-C", repo] + list(a), capture_output=True,
                                  text=True, check=True).stdout
        except (OSError, subprocess.CalledProcessError):
            return None
    head = g("rev-parse", "HEAD")
    status = g("status", "--porcelain", "--untracked-files=all")
    return {"git_commit_submit": os.environ.get("GIT_COMMIT") or None,
            "git_commit_runtime": head.strip() if head else "unavailable",
            "git_dirty": (bool(status.strip()) if status is not None else "unavailable")}


# ------------------------------------------------------------- graph building
def generator(tool, exe, version, tool_args):
    """What makes a cached file (B-1, G-4): the tool, its version as measured
    in this run, the sha256 of the executable that runs, and every argument
    except the file operands. --reuse requires an exact match, so a cache made
    by another build, version or parameter set is never reused."""
    path = shutil.which(exe)
    if not path:
        raise ArmBFail("cannot find executable %r" % exe)
    return {"tool": tool, "version": version, "executable": file_identity(path)[0],
            "executable_sha256": file_identity(path)[2], "args": list(tool_args)}


def creation_record():
    """Who made a cached file; carried into the marker when it is reused."""
    rec = {"utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
           "git_commit_submit": os.environ.get("GIT_COMMIT")}
    rec.update(run_context())
    return rec


def same_generator(meta, gen):
    """Executable path may differ (a moved binary); content and args may not."""
    old = (meta or {}).get("generator") or {}
    keys = ("tool", "version", "executable_sha256", "args")
    return all(k in old for k in keys) and all(old[k] == gen[k] for k in keys)


def build_graph(paths, gfa, mg_args, minigraph):
    """-cxggs = incremental graph construction; the FIRST fasta is the backbone.
    Written to a temporary name and moved into place only on success."""
    tmp = gfa + ".tmp"
    with open(tmp, "w") as fh:
        run([minigraph] + mg_args + paths, stdout=fh, capture_output=False)
    os.replace(tmp, gfa)


def load_segments(gfa):
    segs = {}
    with open(gfa) as fh:
        for line in fh:
            if line.startswith("S\t"):
                f = line.rstrip("\n").split("\t")
                if len(f) < 3 or not f[1] or f[2] in ("", "*"):
                    raise ArmBFail("graph %s has an S line without a sequence: %r"
                                   % (gfa, line[:80]))
                segs[f[1]] = f[2].upper()
    if not segs:
        raise ArmBFail("%s contains no S lines -- minigraph produced no graph" % gfa)
    return segs


def node_seq(token, segs):
    """'>s12' -> segment s12 forward, '<s12' -> reverse complement."""
    m = WALK_TOKEN.fullmatch(token or "")
    if not m or m.group(2) not in segs:
        raise ArmBFail("unknown graph node %r" % token)
    s = segs[m.group(2)]
    return s if m.group(1) == ">" else revcomp(s)


def walk_seq(walk, segs):
    """Sequence of a --call walk ('>s2>s3', '<s4', or '*' for no segment)."""
    if walk == "*":
        return ""
    toks = WALK_TOKEN.findall(walk)
    if not toks or "".join(o + n for o, n in toks) != walk:
        raise ArmBFail("unparseable walk %r" % walk)
    return "".join(node_seq(o + n, segs) for o, n in toks)


def bubbles(gfa, raw, gfatools):
    """`gfatools bubble` emits 14 columns:

        0 chrom  1 start  2 end  3 n_seg  4 n_path  5 is_inversion
        6 len_shortest  7 len_longest  8-10 (-1 unless -d)  11 segment list
        12 shortest-path sequence      13 longest-path sequence

    All 14 columns are checked: 1-10 are integers (8-10 are -1 unless -d),
    11 is a comma-separated segment list, 12/13 are the path sequences, '*'
    for a zero-length path and otherwise exactly len_shortest / len_longest
    letters. Any line that does not match is a parse failure (B-2): returns
    the rows and the offending line numbers, and the caller stops the run."""
    tmp = raw + ".tmp"
    with open(tmp, "w") as fh:
        run([gfatools, "bubble", gfa], stdout=fh, capture_output=False)
    os.replace(tmp, raw)
    rows, bad = {}, []
    with open(raw) as fh:
        for n, line in enumerate(fh, 1):
            f = line.rstrip("\n").split("\t")
            if len(f) != 14:
                bad.append(n)
                continue
            try:
                vals = [int(x) for x in f[1:11]]
            except ValueError:
                bad.append(n)
                continue
            st, en, n_seg, n_path, inv, lo, hi = vals[:7]

            def seq_ok(x, ln):
                # a zero-length path is written '*'; otherwise the sequence
                # must be letters and exactly as long as its length column
                return x == "*" if ln == 0 else bool(SEQ.fullmatch(x)) and len(x) == ln
            if (not f[0] or st < 0 or en < st or n_seg < 1 or n_path < 1
                    or inv not in (0, 1) or lo < 0 or hi < lo
                    or not SEGLIST.fullmatch(f[11])
                    or not seq_ok(f[12], lo) or not seq_ok(f[13], hi)
                    or (f[0], st, en) in rows):
                bad.append(n)
                continue
            rows[(f[0], st, en)] = {"bid": "b%06d" % (n - 1), "n_seg": n_seg,
                                    "n_path": n_path, "inv": inv, "lo": lo, "hi": hi,
                                    "segs": f[11]}
    return rows, bad


# --------------------------------------------------------------- genotyping
def call_sample(gfa, fasta, out_bed, mg_args, minigraph):
    tmp = out_bed + ".tmp"
    with open(tmp, "w") as fh:
        run([minigraph] + mg_args + [gfa, fasta], stdout=fh, capture_output=False)
    os.replace(tmp, out_bed)


MISSING_CALL = (None, None, None, None, ".", None, None, None)


def parse_call(bed):
    """`minigraph -cxasm --call` emits 6 columns:

        0 chrom  1 start  2 end  3 source node  4 sink node
        5 walk : allele_len : strand : query_contig : query_start : query_end

    e.g.  >s2:128:-:CP154596.1:4269163:4269298
    A missing call ('*', '.') is MISSING_CALL; a zero-length allele has
    allele_len 0 and keeps its query interval. Anything else is a parse
    failure (B-2): no field is ever guessed.
    Returns ({(chrom, start, end): (allele_len, contig, qs, qe, strand, walk,
    source, sink)}, [bad line numbers]).
    """
    out, bad = {}, []
    with open(bed) as fh:
        for n, line in enumerate(fh, 1):
            f = line.rstrip("\n").split("\t")
            if len(f) != 6:
                bad.append(n)
                continue
            try:
                key = (f[0], int(f[1]), int(f[2]))
            except ValueError:
                bad.append(n)
                continue
            if key in out:
                bad.append(n)                 # a bubble may appear once per file
                continue
            if not (WALK_TOKEN.fullmatch(f[3]) and WALK_TOKEN.fullmatch(f[4])):
                bad.append(n)                 # source / sink must be one node each
                continue
            rec = f[5].strip()
            if rec in MISSING:
                out[key] = MISSING_CALL
                continue
            tok = rec.split(":")
            if len(tok) != 6:
                bad.append(n)
                continue
            try:
                alen, qs, qe = int(tok[1]), int(tok[4]), int(tok[5])
            except ValueError:
                bad.append(n)
                continue
            walk_ok = tok[0] == "*" or "".join(
                o + x for o, x in WALK_TOKEN.findall(tok[0])) == tok[0]
            if (tok[2] not in ("+", "-") or not tok[3] or not tok[0] or not walk_ok
                    or (tok[0] == "*") != (alen == 0)
                    or alen < 0 or qs < 0 or qe < qs):
                bad.append(n)
                continue
            out[key] = (alen, tok[3], qs, qe, tok[2], tok[0], f[3], f[4])
    return out, bad


# ------------------------------------------------------------- classification
def classify(obs, len_tol, min_insert):
    """THE ONLY place the allele classes are defined (B-3).

    obs: {sample: allele_len}. empty = within len_tol of the shortest allele;
    filled = at least min_insert longer than the shortest; intermediate = in
    between. Requires len_tol < min_insert (checked at start), so the classes
    are disjoint and exhaustive over the observed genomes."""
    shortest = min(obs.values())
    empty = sorted(s for s, a in obs.items() if a <= shortest + len_tol)
    filled = sorted(s for s, a in obs.items() if a >= shortest + min_insert)
    inter = sorted(s for s in obs if s not in empty and s not in filled)
    return shortest, empty, filled, inter


def length_classes(lengths, len_tol):
    """Number of allele-length classes. Each class is anchored at its SHORTEST
    member and holds lengths within len_tol of that anchor, so intermediate
    lengths cannot chain classes together (0, 40, 80, 800 -> 3 classes)."""
    n, anchor = 0, None
    for x in sorted(lengths):
        if anchor is None or x - anchor > len_tol:
            n, anchor = n + 1, x
    return n


def find_hits(hay, needle, cap=2):
    """Start positions of needle in hay (overlapping hits counted), at most cap."""
    hits, i = [], hay.find(needle)
    while i != -1 and len(hits) < cap:
        hits.append(i)
        i = hay.find(needle, i + 1)
    return hits


def locate_allele(contig_seq, call, segs, context):
    """B-6: place a genome's OWN allele in its contig.

    The needle is source-node tail + the genome's walk + sink-node head, in
    graph orientation. It is accepted as the frame only if
      (1) it occurs exactly once in the WHOLE contig, counting both strands,
      (2) that occurrence is on the call's strand, and
      (3) it overlaps the minigraph call interval of this genome (by at least
          half of the shorter of the two spans; a zero-length call interval
          must fall inside it), and
      (4) the genome's whole call interval equals the graph path (source
          node + walk + sink node) at the offset this occurrence implies,
    so a repeat copy -- outside the interval or inside it -- can never stand
    in for the called site.
    Returns (status, a0, a1): the allele interior as 0-based half-open contig
    coordinates, or (reason, None, None) with reason in unverified (no exact
    occurrence), ambiguous (several, or on the other strand),
    off_call_interval (one exact occurrence, but not at the called site) or
    interval_mismatch (the call interval's own sequence does not line up)."""
    alen, contig, qs, qe, strand, walk, source, sink = call
    ws = walk_seq(walk, segs)
    left, right = node_seq(source, segs)[-context:], node_seq(sink, segs)[:context]
    needle = left + ws + right
    probe = needle if strand == "+" else revcomp(needle)
    other = revcomp(probe)
    hits = find_hits(contig_seq, probe)
    other_hits = [] if other == probe else find_hits(contig_seq, other, cap=1)
    if not hits and not other_hits:
        return "unverified", None, None
    if len(hits) != 1 or other_hits:
        return "ambiguous", None, None
    n0, n1 = hits[0], hits[0] + len(needle)
    if qe > qs:
        ov = max(0, min(n1, qe) - max(n0, qs))
        at_call = ov >= 0.5 * min(n1 - n0, qe - qs)
    else:
        at_call = n0 <= qs <= n1
    if not at_call:
        return "off_call_interval", None, None
    a0 = n0 + (len(left) if strand == "+" else len(right))
    a1 = a0 + len(ws)
    # (4) the genome's WHOLE call interval must equal the graph path
    # source + walk + sink at the offset this occurrence implies. A copy that
    # merely lies inside a wide interval while the called site differs (say by
    # a SNP) cannot pass: the interval content would not line up.
    src, snk = node_seq(source, segs), node_seq(sink, segs)
    path = src + ws + snk
    iv = contig_seq[qs:qe]
    if strand == "+":
        ps = len(src) - (a0 - qs)
    else:
        iv = revcomp(iv)
        ps = len(src) - (qe - a1)
    if ps < 0 or ps + len(iv) > len(path) or path[ps:ps + len(iv)] != iv:
        return "interval_mismatch", None, None
    return "verified", a0, a1


def oriented(fa, contig, start, end, strand):
    """Contig slice [start, end) in graph orientation."""
    s = fa.fetch(contig, max(0, start), min(fa.length(contig), end)).upper()
    return s if strand == "+" else revcomp(s)


# ------------------------------------------------------------- main
def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", required=True,
                    help="panel manifest, one FASTA path per line (backbone on line 1)")
    ap.add_argument("--min-genome-frac", type=float, default=0.5,
                    help="reject a manifest entry smaller than this fraction of the "
                         "panel's MEDIAN genome size (plasmid-only 'Complete Genome' "
                         "entries; relative, so small genomes are not rejected)")
    ap.add_argument("--min-genome-bytes", type=int, default=100000,
                    help="absolute floor, for a file that is not an assembly at all")
    ap.add_argument("--allow-zero-bubbles", action="store_true",
                    help="record a genuine zero-bubble panel instead of failing")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--cluster-id", default=None)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--min-insert", type=int, default=500,
                    help="min allele-length difference to call a structural allele")
    ap.add_argument("--max-insert", type=int, default=5000,
                    help="primary scope; larger events are parked, not discarded")
    ap.add_argument("--len-tol", type=int, default=50,
                    help="allele lengths within this of the shortest are empty; "
                         "must be < --min-insert")
    ap.add_argument("--min-empty", type=int, default=1)
    ap.add_argument("--min-filled", type=int, default=1)
    ap.add_argument("--flank-qc", type=int, default=500,
                    help="required unique anchor each side, inside one contig; "
                         "also the flank exported around empty sites")
    ap.add_argument("--clean-span-tol", type=int, default=20,
                    help="junction spans within this are reported as "
                         "clean_or_graph_slack, not asymmetric")
    ap.add_argument("--locate-context", type=int, default=50,
                    help="bp of source/sink node sequence around a walk when "
                         "placing it in its genome (B-6)")
    ap.add_argument("--min-tsd", type=int, default=3)
    ap.add_argument("--max-tsd", type=int, default=30)
    ap.add_argument("--keep-inversions", action="store_true")
    ap.add_argument("--minigraph", default="minigraph")
    ap.add_argument("--gfatools", default="gfatools")
    ap.add_argument("--samtools", default="samtools",
                    help="used only to build a missing .fai; its version is then recorded")
    ap.add_argument("--reuse", action="store_true",
                    help="reuse graph/calls, but only when bound to the same inputs")
    return ap.parse_args(argv)


def write_failure(outdir, code, reason, details):
    """B-2 / B-1: a failed run leaves a failure summary (never a marker)."""
    try:
        rec = {"status": "failed", "exit_code": code, "reason": reason,
               "details": details, "argv": sys.argv}
        rec.update(run_context())
        rec.update(git_provenance())
        write_json(os.path.join(outdir, FAILED), rec)
    except Exception:
        pass        # the FATAL log line and exit code remain the primary signal


def main(argv=None):
    args = parse_args(argv)
    tmpdir = None
    lock_fh = None
    locked = False
    try:
        # ---- B-1: one run per output directory, and the old completion
        # marker is invalidated before ANY check that can fail -------------
        os.makedirs(args.outdir, exist_ok=True)
        lock_fh = open(os.path.join(args.outdir, LOCKFILE), "w")
        try:
            fcntl.flock(lock_fh, fcntl.LOCK_EX | fcntl.LOCK_NB)
            locked = True
        except OSError:
            raise ArmBFail("another Arm B run holds %s; one run per output directory"
                           % os.path.join(args.outdir, LOCKFILE))
        for name in (MARKER, FAILED) + OUTPUTS + ("bubbles.tsv.raw",):
            p = os.path.join(args.outdir, name)
            if os.path.lexists(p):
                os.remove(p)
        for d in os.listdir(args.outdir):      # leftovers of crashed runs
            if d.startswith(TMPDIR_PREFIX):
                shutil.rmtree(os.path.join(args.outdir, d), ignore_errors=True)
        tmpdir = tempfile.mkdtemp(prefix=TMPDIR_PREFIX, dir=args.outdir)
        return run_armb(args, tmpdir)
    except ArmBFail as e:
        log("FATAL: %s" % e)
        if tmpdir:
            shutil.rmtree(tmpdir, ignore_errors=True)
        if locked:      # only the lock holder may touch the shared directory
            write_failure(args.outdir, e.code, str(e), e.details)
        return e.code
    except subprocess.CalledProcessError as e:
        log("FATAL: external tool failed (exit %s): %s"
            % (e.returncode, " ".join(map(str, e.cmd))))
        if tmpdir:
            shutil.rmtree(tmpdir, ignore_errors=True)
        if locked:
            write_failure(args.outdir, 2, "external tool failed (exit %s)" % e.returncode,
                          {"cmd": list(map(str, e.cmd))})
        return 2
    except Exception as e:      # OSError, short reads, missing executables, bugs
        import traceback
        traceback.print_exc()
        log("FATAL: %s: %s" % (type(e).__name__, e))
        if tmpdir:
            shutil.rmtree(tmpdir, ignore_errors=True)
        if locked:
            write_failure(args.outdir, 2, "%s: %s" % (type(e).__name__, e), {})
        return 2
    finally:
        if lock_fh:
            lock_fh.close()          # releases the flock


def run_armb(args, tmpdir):
    # ---- parameters ---------------------------------------------------------
    if not 0 <= args.len_tol < args.min_insert <= args.max_insert:
        raise ArmBFail("need 0 <= --len-tol < --min-insert <= --max-insert (got %d, %d, %d)"
                       % (args.len_tol, args.min_insert, args.max_insert))
    if args.locate_context < 1 or args.flank_qc < 0:
        raise ArmBFail("--locate-context must be >= 1 and --flank-qc >= 0")
    if not re.fullmatch(r"[0-9a-f]{40}", os.environ.get("GIT_COMMIT", "")):
        raise ArmBFail("GIT_COMMIT is not set to a 40-hex commit (G-4): run as "
                       "GIT_COMMIT=$(git rev-parse HEAD) ... or sbatch --export=ALL,GIT_COMMIT=...")
    versions = tool_versions(args.minigraph, args.gfatools)
    cid = args.cluster_id or os.path.basename(args.manifest).split(".")[0]
    os.makedirs(os.path.join(args.outdir, "calls"), exist_ok=True)

    # ---- inputs: exist, are assemblies, and are distinct --------------------
    paths = [l.strip() for l in open(args.manifest) if l.strip()]
    if len(paths) < 2:
        raise ArmBFail("manifest lists %d genome(s); a panel needs at least 2" % len(paths))
    missing = [p for p in paths if not os.path.exists(p)]
    if missing:
        raise ArmBFail("%d of %d manifest entries do not exist, first: %s"
                       % (len(missing), len(paths), missing[0]))
    real = [os.path.realpath(p) for p in paths]
    names = [sample_name(p) for p in paths]
    dup_paths = sorted({p for p in real if real.count(p) > 1})
    dup_names = sorted({n for n in names if names.count(n) > 1})
    if dup_paths or dup_names:
        raise ArmBFail("manifest entries must be distinct files with distinct sample names; "
                       "duplicate files %s, duplicate names %s" % (dup_paths[:3], dup_names[:3]))
    sizes = [os.path.getsize(p) for p in paths]
    med = sorted(sizes)[len(sizes) // 2]
    floor = max(args.min_genome_bytes, int(args.min_genome_frac * med))
    tiny = [(p, z) for p, z in zip(paths, sizes) if z < floor]
    if tiny:
        raise ArmBFail("%d of %d manifest entries are below %d bytes (%.0f%% of the panel "
                       "median %d); not whole assemblies. First: %s (%d bytes)"
                       % (len(tiny), len(paths), floor, 100 * args.min_genome_frac, med,
                          tiny[0][0], tiny[0][1]))
    log("%s: %d genomes, backbone=%s; size floor %d bytes, all pass"
        % (cid, len(paths), names[0], floor))
    need_index = [p for p in paths if not os.path.exists(p + ".fai")]
    if need_index:
        out = run([args.samtools, "--version"]).stdout.strip().splitlines()
        if not out:
            raise ArmBFail("samtools printed no version")
        versions["samtools"] = out[0]
        log("building %d missing .fai index(es) with %s" % (len(need_index), out[0]))
    else:
        versions["samtools"] = "not invoked (every .fai already present)"
    fais = {}
    for p, sn in zip(paths, names):
        try:
            fais[sn] = Fasta(p, samtools=args.samtools)
        except Exception as e:      # index build/read failure is fatal (B-1)
            raise ArmBFail("cannot index/read %s: %s" % (p, e))

    # ---- graph, bound to its inputs and its generator for --reuse ------------
    gfa = os.path.join(args.outdir, "graph.gfa")
    gmeta_path = gfa + ".meta.json"
    inputs_id = [file_identity(p) for p in paths]
    graph_args = ["-cxggs", "-t", str(args.threads)]
    call_args = ["-cxasm", "--call", "-t", str(args.threads)]
    graph_gen = generator("minigraph", args.minigraph, versions["minigraph"], graph_args)
    call_gen = generator("minigraph", args.minigraph, versions["minigraph"], call_args)
    gmeta = read_json(gmeta_path) if args.reuse else None
    reuse_graph = bool(gmeta and gmeta.get("inputs") == inputs_id and os.path.exists(gfa)
                       and same_generator(gmeta, graph_gen)
                       and gmeta.get("gfa_sha256") == sha256_file(gfa))
    if reuse_graph:
        log("reusing graph (same inputs and generator, sha256 verified)")
    else:
        if os.path.lexists(gmeta_path):
            os.remove(gmeta_path)
        log("building graph ...")
        build_graph(paths, gfa, graph_args, args.minigraph)
        if not os.path.exists(gfa) or os.path.getsize(gfa) == 0:
            raise ArmBFail("%s is missing or empty after graph construction" % gfa)
        gmeta = {"inputs": inputs_id, "gfa_sha256": sha256_file(gfa),
                 "generator": graph_gen, "created": creation_record()}
        write_json(gmeta_path, gmeta)
    gfa_sha = sha256_file(gfa)
    cache = {"graph": {"reused": reuse_graph, "meta": gmeta}, "calls": {}}
    segs = load_segments(gfa)
    log("graph: %.1f MB, %d segments" % (os.path.getsize(gfa) / 1e6, len(segs)))

    bindex, bad = bubbles(gfa, os.path.join(tmpdir, "bubbles.raw"), args.gfatools)
    if bad:
        raise ArmBFail("%d unparseable gfatools bubble line(s) (expected 14 documented "
                       "columns), first at line %d" % (len(bad), bad[0]),
                       details={"unparseable_bubble_lines": len(bad), "first_bad_line": bad[0]})
    with open(os.path.join(tmpdir, "bubbles.tsv"), "w") as out:
        out.write("bubble_id\tchrom\tstart\tend\tn_seg\tn_path\tis_inversion\t"
                  "len_shortest\tlen_longest\tlen_delta\tsegments\n")
        for (chrom, st, en), b in bindex.items():
            out.write("%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n"
                      % (b["bid"], chrom, st, en, b["n_seg"], b["n_path"], b["inv"],
                         b["lo"], b["hi"], b["hi"] - b["lo"], b["segs"]))
    log("%d bubbles" % len(bindex))

    # ---- genotype every genome, call files bound to graph + FASTA ------------
    calls, walk_len = {}, {}
    for p, sn in zip(paths, names):
        bed = os.path.join(args.outdir, "calls", sn + ".bed")
        bmeta_path = bed + ".meta.json"
        bmeta = read_json(bmeta_path) if args.reuse else None
        ok = bool(bmeta and bmeta.get("gfa_sha256") == gfa_sha
                  and bmeta.get("fasta") == file_identity(p) and os.path.exists(bed)
                  and same_generator(bmeta, call_gen)
                  and bmeta.get("bed_sha256") == sha256_file(bed))
        if not ok:
            if os.path.lexists(bmeta_path):
                os.remove(bmeta_path)
            call_sample(gfa, p, bed, call_args, args.minigraph)
            bmeta = {"gfa_sha256": gfa_sha, "fasta": file_identity(p),
                     "bed_sha256": sha256_file(bed), "generator": call_gen,
                     "created": creation_record()}
            write_json(bmeta_path, bmeta)
        cache["calls"][sn] = {"reused": ok, "meta": bmeta}
        calls[sn], bad = parse_call(bed)
        if bad:
            raise ArmBFail("%d unparseable minigraph --call line(s) in %s (expected 6 "
                           "documented fields), first at line %d" % (len(bad), bed, bad[0]),
                           details={"unparseable_call_lines": len(bad), "file": bed,
                                    "first_bad_line": bad[0]})
        # minigraph --call writes one line per bubble, '*' when the genome has
        # no allele there. An absent RECORD is therefore a truncated or wrong
        # file, never a missing call.
        # B-1/B-2: every non-missing call is checked now, before any filtering
        # or choice of reference genome: contig present, interval inside it,
        # source/sink/walk nodes in the graph, walk length == allele length
        fa = fais[sn]
        for k, c in calls[sn].items():
            if c[0] is None:
                continue
            alen, contig, qs, qe, strand, walk, source, sink = c
            if contig not in fa:
                raise ArmBFail("%s: call contig %r at %s:%d-%d is not in %s's FASTA"
                               % (bed, contig, k[0], k[1], k[2], sn))
            if qe > fa.length(contig):
                raise ArmBFail("%s: call interval %d-%d runs past the end of %s (%d bp)"
                               % (bed, qs, qe, contig, fa.length(contig)))
            node_seq(source, segs)
            node_seq(sink, segs)
            if walk not in walk_len:
                walk_len[walk] = len(walk_seq(walk, segs))
            if walk_len[walk] != alen:
                raise ArmBFail("%s: walk %s has %d bp but the call says %d at %s:%d-%d"
                               % (bed, walk[:60], walk_len[walk], alen, k[0], k[1], k[2]))
        absent = [k for k in bindex if k not in calls[sn]]
        extra = [k for k in calls[sn] if k not in bindex]
        if absent or extra:
            raise ArmBFail("%s does not cover the graph's bubbles: %d bubble(s) without a "
                           "record, %d record(s) for unknown bubbles (first: %s)"
                           % (bed, len(absent), len(extra), (absent or extra)[0]))
        log("  called %s%s (%d bubbles covered)"
            % (sn, " [reused]" if ok else "",
               sum(1 for v in calls[sn].values() if v[0] is not None)))

    # ---- classify bubbles ----------------------------------------------------
    T = {k: open(os.path.join(tmpdir, k), "w") for k in OUTPUTS if k != "bubbles.tsv"}
    T["armB_events.tsv"].write("\t".join(EVENT_COLS) + "\n")
    T["armB_large_events.tsv"].write("\t".join(LARGE_COLS) + "\n")
    T["armB_pairs.tsv"].write("event_id\tcluster_id\tfilled_genome\tempty_genome\tchrom\t"
                              "bubble_start\tbubble_end\tinsert_size\n")

    contig_cache = {}

    def contig_seq(sn, contig):
        k = (sn, contig)
        if k not in contig_cache:
            contig_cache[k] = fais[sn].fetch(contig).upper()
        return contig_cache[k]

    def contiguity(sn):
        fa = fais[sn]
        return -(fa.total() / max(1, len(fa.names())))

    counts = dict(bubbles=len(bindex), unparseable_bubble_lines=0, unparseable_call_lines=0,
                  skipped_inversion=0, fewer_than_2_called=0,
                  too_few_empty_or_filled=0, below_min_insert=0, primary=0, parked=0,
                  frame_verified=0, frame_unverified=0, frame_ambiguous=0,
                  frame_off_call_interval=0, frame_interval_mismatch=0,
                  insert_len_mismatch=0, multiallelic=0)
    n_ev = n_large = 0
    for key, b in bindex.items():
        chrom, bs, be = key
        if b["inv"] and not args.keep_inversions:
            counts["skipped_inversion"] += 1
            continue
        obs, n_missing = {}, 0
        for sn in names:
            c = calls[sn][key]
            if c[0] is None:
                n_missing += 1
            else:
                obs[sn] = c
        if len(obs) < 2:
            counts["fewer_than_2_called"] += 1
            continue
        alens = {s: c[0] for s, c in obs.items()}
        shortest, empty, filled, inter = classify(alens, args.len_tol, args.min_insert)
        if not filled:
            counts["below_min_insert"] += 1         # no allele reaches the size gate
            continue
        if len(empty) < args.min_empty or len(filled) < args.min_filled:
            counts["too_few_empty_or_filled"] += 1
            continue
        len_filled = statistics.median(alens[s] for s in filled)
        insert = int(len_filled - statistics.median(alens[s] for s in empty))
        if insert < args.min_insert:
            counts["below_min_insert"] += 1
            continue
        n_cls = length_classes(alens.values(), args.len_tol)
        flags = []
        if n_missing:
            flags.append("missing_in_%d" % n_missing)
        if inter:
            flags.append("intermediate_%d" % len(inter))
        if n_cls > 2 or inter:
            flags.append("multiallelic")
            counts["multiallelic"] += 1
        if b["inv"]:
            flags.append("inversion")

        # B-5: the two alleles are the called genomes' own walks
        frep = min(filled, key=lambda s: (abs(alens[s] - len_filled), contiguity(s), s))
        ordered_empty = sorted(empty, key=lambda s: (contiguity(s), s))
        frame, eref, a0, a1 = None, ordered_empty[0], None, None
        if insert <= args.max_insert:
            reasons = []
            for s in ordered_empty:            # prefer an empty genome that verifies
                st, x0, x1 = locate_allele(contig_seq(s, obs[s][1]), obs[s], segs,
                                           args.locate_context)
                reasons.append(st)
                if st == "verified":
                    frame, eref, a0, a1 = st, s, x0, x1
                    break
            if frame is None:
                # no empty genome verified: report the declared reference (the
                # most contiguous) with ITS OWN reason, never another genome's
                eref, frame = ordered_empty[0], reasons[0]
        empty_ws = walk_seq(obs[eref][5], segs)
        filled_ws = walk_seq(obs[frep][5], segs)
        if len(empty_ws) != alens[eref] or len(filled_ws) != alens[frep]:
            raise ArmBFail("walk length disagrees with called allele length at %s:%d-%d"
                           % (chrom, bs, be))
        ins_seq, aL, aR = diff_alleles(empty_ws, filled_ws)
        if abs(len(ins_seq) - insert) > args.len_tol:
            flags.append("insert_len_mismatch")
            counts["insert_len_mismatch"] += 1

        if insert > args.max_insert:
            # parked: out of the primary catalogue, kept with its sequence
            lid = "%s.L%06d" % (cid, n_large)
            T["armB_large_events.tsv"].write("\t".join(map(str, [
                lid, cid, chrom, bs, be, b["inv"], shortest, int(len_filled), insert,
                len(empty), len(filled), len(inter), n_missing, n_cls,
                eref, frep, len(ins_seq), seq_md5(ins_seq) if ins_seq else "NA",
                ",".join(empty[:20]), ",".join(filled[:20]), ";".join(flags) or ".",
            ])) + "\n")
            if ins_seq:
                T["armB_large_inserts.fna"].write(">%s len=%d\n%s\n" % (lid, len(ins_seq), ins_seq))
            n_large += 1
            counts["parked"] += 1
            continue

        # FNA-specific QC: does a single contig span the locus in every genome?
        tiers = []
        for s in filled + empty:
            alen, contig, qs, qe = obs[s][:4]
            clen = fais[s].length(contig)
            left, right = qs, clen - qe
            tiers.append("high" if left >= args.flank_qc and right >= args.flank_qc
                         else "medium" if left >= 100 and right >= 100
                         else "unresolvable_contig_break")
        tier = ("high" if all(t == "high" for t in tiers)
                else "unresolvable_contig_break" if "unresolvable_contig_break" in tiers
                else "medium")

        e_call = obs[eref]
        e_contig, e_qs, e_qe, e_strand = e_call[1], e_call[2], e_call[3], e_call[4]
        efa = fais[eref]
        if frame == "verified":
            counts["frame_verified"] += 1
            # B-6: interior [a0, a1) of eref's own allele; anchors from its walk
            if e_strand == "+":
                jL, jR = a0 + aL, a1 - aR
            else:
                jL, jR = a0 + aR, a1 - aL
            span = jR - jL
            arch = "clean_or_graph_slack" if abs(span) <= args.clean_span_tol \
                else "asymmetric_candidate"
            # flanks from the genome itself, in graph orientation
            if e_strand == "+":
                lctx = oriented(efa, e_contig, jL - args.flank_qc, jL, "+")
                rctx = oriented(efa, e_contig, jR, jR + args.flank_qc, "+")
                site = oriented(efa, e_contig, jL - args.flank_qc, jR + args.flank_qc, "+")
            else:
                lctx = oriented(efa, e_contig, jR, jR + args.flank_qc, "-")
                rctx = oriented(efa, e_contig, jL - args.flank_qc, jL, "-")
                site = oriented(efa, e_contig, jL - args.flank_qc, jR + args.flank_qc, "-")
        else:
            counts["frame_" + frame] += 1
            jL = jR = span = -1
            arch = "no_empty_frame"
            flags.append("frame_" + frame)
            src = node_seq(e_call[6], segs)[-args.flank_qc:]
            snk = node_seq(e_call[7], segs)[:args.flank_qc]
            lctx = src + empty_ws[:aL]
            rctx = empty_ws[len(empty_ws) - aR:] + snk
            site = src + empty_ws + snk
        # B-7: TSD is recorded, never used to decide anything
        tl, ts, tside, tconf = find_tsd(ins_seq, lctx, rctx,
                                       min_tsd=args.min_tsd, max_tsd=args.max_tsd)
        eid = "%s.B%06d" % (cid, n_ev)
        T["armB_events.tsv"].write("\t".join(map(str, [
            eid, cid, chrom, bs, be, b["inv"], shortest, int(len_filled), insert,
            len(empty), len(filled), len(inter), n_missing, n_cls,
            "PRESENCE_ABSENCE", eref, e_contig, e_strand, e_qs, e_qe, frame,
            jL, jR, span, arch, frep, len(ins_seq),
            seq_md5(ins_seq) if ins_seq else "NA",
            tl, ts or ".", tside, tconf,
            ",".join(empty[:20]), ",".join(filled[:20]), tier, ";".join(flags) or ".",
        ])) + "\n")
        if ins_seq:
            T["armB_inserts.fna"].write(">%s len=%d tsd=%d\n%s\n" % (eid, len(ins_seq), tl, ins_seq))
        T["armB_empty_sites.fna"].write(">%s len=%d frame=%s\n%s\n" % (eid, len(site), frame, site))
        for f in sorted(filled, key=lambda s: (contiguity(s), s))[:max(1, args.min_filled)]:
            for e in ordered_empty[:max(1, args.min_empty)]:
                T["armB_pairs.tsv"].write("%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\n"
                                          % (eid, cid, f, e, chrom, bs, be, insert))
        n_ev += 1
        counts["primary"] += 1
    for fh in T.values():
        fh.close()

    log("%s: %s" % (cid, json.dumps(counts, sort_keys=True)))

    # ---- B-9: zero bubbles is a failure until proven otherwise ---------------
    if len(bindex) == 0 and not args.allow_zero_bubbles:
        raise ArmBFail("%s produced 0 bubbles from %d genomes. Either the panel is "
                       "genuinely invariant or the run failed silently -- these look the "
                       "same. Re-run with --allow-zero-bubbles to record a genuine zero."
                       % (cid, len(paths)), code=3)

    # ---- publish: move every output into place, then the marker last ---------
    shas = {}
    for name in OUTPUTS:
        src = os.path.join(tmpdir, name)
        shas[name] = sha256_file(src)
        os.replace(src, os.path.join(args.outdir, name))
    os.replace(os.path.join(tmpdir, "bubbles.raw"), os.path.join(args.outdir, "bubbles.tsv.raw"))
    shutil.rmtree(tmpdir, ignore_errors=True)
    marker = {"status": "complete", "cluster_id": cid, "argv": sys.argv,
              "params": {k: v for k, v in vars(args).items()},
              "tool_versions": versions, "graph_sha256": gfa_sha,
              # B-1/G-4: for the graph and each call file, whether it was
              # reused and the generator + creation record of the run that
              # made it (a reused file keeps its original record)
              "cache": cache,
              "manifest": {"path": os.path.realpath(args.manifest),
                           "sha256": sha256_file(args.manifest)},
              "inputs": [{"sample": sn, "path": file_identity(p)[0],
                          "bytes": file_identity(p)[1], "sha256": file_identity(p)[2]}
                         for p, sn in zip(paths, names)],
              "counts": counts, "outputs_sha256": shas}
    marker.update(run_context())
    marker.update(git_provenance())
    write_json(os.path.join(args.outdir, MARKER), marker)
    log("%s: %d empty/filled structural alleles, %d parked (>%d bp); complete"
        % (cid, n_ev, n_large, args.max_insert))
    return 0


if __name__ == "__main__":
    sys.exit(main())
