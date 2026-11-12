#!/usr/bin/env python3
"""M2 corpus: rebuild the corpus table, map duplicate assemblies, report gaps.

SPEC item M2-C (docs/SPEC.md). Three subcommands. Each writes into a NEW
output directory (an existing one is refused); outputs are staged in
<out>/.tmp and moved into place only after the whole run succeeded, and the
run record <cmd>.json (G-4: commits, host, SLURM ids, argv, input sha256,
counts, output sha256) is written LAST. A failed run leaves FAILED.json and
no outputs.

  build   config/corpus_sources.tsv -> corpus.tsv
          Every *.fna at the top level of each declared source directory
          (sub-directories such as _plasmid_only/ are not entered), labelled
          with the directory's species, sorted by accession. Each directory
          must hold exactly its declared expected_n files. --expect-sha256
          pins the result (corpus_v2.tsv must be rebuilt byte for byte).

  pairs   corpus.tsv + GenBank summary
          -> corpus_content.tsv, corpus_duplicates.tsv, corpus_identical_other_id.tsv
          Reads EVERY genome of the corpus once, completely (bytes read must
          equal the file size), and records each file's path, bytes, sha256
          and sequence content in corpus_content.tsv -- the full list of the
          inputs this map was made from.
          (1) Assemblies sharing one 9-digit accession id (GCA_x / GCF_x, or
          two versions) are compared by sequence content. Members that are
          the same assembly (identical sequence, or a different version of
          it) are marked collapse=yes into one canonical accession: counts
          per isolate use this map (user decision 2026-11-12, 1(a)).
          (2) Assemblies with DIFFERENT ids but identical sequence content
          are listed only, never collapsed: identical sequence can come from
          distinct isolates of a clonal species.

  gaps    corpus.tsv + config + GenBank summary + cutoff date
          -> corpus_gaps.tsv (per species), corpus_gap_list.tsv
          Usable = assembly_level Complete Genome or Chromosome, matched on
          species_taxid. Matching is by 9-digit accession id.

Sequence content of an assembly (pairs): every record's sequence with all
ASCII whitespace (space, tab, LF, CR, VT, FF) removed, upper-cased, sha256 per record; the assembly's
seq_sha256 is the sha256 of the sorted record digests. Headers, record order,
line width, case and line endings do not matter; one base does.
"""
import argparse
import collections
import csv
import hashlib
import json
import multiprocessing
import os
import re
import shutil
import socket
import stat
import subprocess
import sys
import time

ACC = re.compile(r"(GC[AF])_(\d{9})\.(\d+)")
DATE = re.compile(r"\d{4}-\d{2}-\d{2}")
HEX64 = re.compile(r"[0-9a-f]{64}")
# an http(s) URL with a host name (optional :port) and a path; the same rule
# as tools/download_staging.sh (tests/test_build_corpus.py checks they agree)
FTP = re.compile(r"https?://[A-Za-z0-9](?:[A-Za-z0-9.-]*[A-Za-z0-9])?(?::[0-9]+)?"
                 r"/[A-Za-z0-9._~/:+-]*")
SOURCES_COLS = ["species", "species_taxid", "organism", "source_dir", "expected_n"]
CORPUS_COLS = ["assembly", "species", "fna_path"]
USABLE_LEVELS = ("Complete Genome", "Chromosome")
SUMMARY_NEED = ("assembly_accession", "species_taxid", "assembly_level", "seq_rel_date",
                "ftp_path", "gbrs_paired_asm", "paired_asm_comp")
PLASMID_DIR = "_plasmid_only"
ASCII_WS = b" \t\n\r\x0b\x0c"       # every byte bytes.isspace() accepts
LIST_TRIES = 4


class Fail(Exception):
    pass


def log(msg):
    sys.stderr.write("[250] %s\n" % msg)
    sys.stderr.flush()


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for blk in iter(lambda: fh.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def acc_parts(a):
    m = ACC.fullmatch(a)
    if not m:
        raise Fail("not an NCBI assembly accession: %r" % a)
    return m.group(1), m.group(2), int(m.group(3))


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
    if not head or not re.fullmatch(r"[0-9a-f]{40}", head.strip()) or status is None:
        raise Fail("cannot measure git_commit_runtime / git_dirty in %s (G-4)" % repo)
    return {"git_commit_submit": os.environ.get("GIT_COMMIT"),
            "git_commit_runtime": head.strip(), "git_dirty": bool(status.strip())}


class Input:
    """One input table, read ONCE: the bytes parsed are exactly the bytes
    counted and hashed (no separate checksum pass that a later read could
    disagree with). When lines() is exhausted, the bytes read must equal the
    file size at open and at EOF, and the sha256 must equal the pin if one is
    given; only then is the record (path, bytes, sha256) set for run.input().
    A malformed pin is fatal before anything is read."""

    def __init__(self, path, what, want=None):
        if want is not None and not HEX64.fullmatch(want):
            raise Fail("%s sha256 must be 64 lowercase hex, got %r" % (what, want))
        self.path, self.what, self.want, self.rec = path, what, want, None

    def lines(self):
        h, n = hashlib.sha256(), 0
        try:
            fh = open(self.path, "rb")
        except OSError as e:
            raise Fail("cannot open %s %s: %s" % (self.what, self.path, e))
        with fh:
            size0 = os.fstat(fh.fileno()).st_size
            for raw in fh:
                h.update(raw)
                n += len(raw)
                try:
                    line = raw.decode("utf-8")
                except UnicodeDecodeError as e:
                    raise Fail("%s %s is not UTF-8 text: %s" % (self.what, self.path, e))
                yield line
            size1 = os.fstat(fh.fileno()).st_size
        if not n == size0 == size1:
            raise Fail("%s %s: short read: %d bytes read, file has %d at open, %d at EOF"
                       % (self.what, self.path, n, size0, size1))
        got = h.hexdigest()
        if self.want is not None and got != self.want:
            raise Fail("%s %s has sha256 %s, expected %s" % (self.what, self.path, got, self.want))
        self.rec = {"path": os.path.realpath(self.path), "bytes": n, "sha256": got}


# ------------------------------------------------------------------ readers
def read_sources(inp):
    path = inp.path
    rows = list(csv.reader(inp.lines(), delimiter="\t"))
    if not rows or rows[0] != SOURCES_COLS:
        raise Fail("%s: header must be %s" % (path, "\t".join(SOURCES_COLS)))
    out, seen_sp, seen_dir = [], set(), set()
    for n, r in enumerate(rows[1:], 2):
        if len(r) != len(SOURCES_COLS):
            raise Fail("%s line %d: %d fields, expected %d" % (path, n, len(r), len(SOURCES_COLS)))
        sp, taxid, org, d, exp = r
        if (not re.fullmatch(r"[a-z][a-z0-9_]*", sp) or not taxid.isdigit()
                or not os.path.isabs(d) or d != os.path.normpath(d) or not exp.isdigit()):
            raise Fail("%s line %d: invalid row %r" % (path, n, r))
        if sp in seen_sp or d in seen_dir:
            raise Fail("%s line %d: duplicate species or source_dir" % (path, n))
        seen_sp.add(sp)
        seen_dir.add(d)
        out.append({"species": sp, "species_taxid": taxid, "organism": org,
                    "source_dir": d, "expected_n": int(exp)})
    if not out:
        raise Fail("%s declares no sources" % path)
    return out


def read_corpus(inp):
    path = inp.path
    rows = list(csv.reader(inp.lines(), delimiter="\t"))
    if not rows or rows[0] != CORPUS_COLS:
        raise Fail("%s: header must be %s" % (path, "\t".join(CORPUS_COLS)))
    out, seen = [], set()
    for n, r in enumerate(rows[1:], 2):
        if len(r) != 3 or not all(r):
            raise Fail("%s line %d: expected 3 non-empty fields" % (path, n))
        acc_parts(r[0])
        if r[0] in seen:
            raise Fail("%s line %d: duplicate assembly %s" % (path, n, r[0]))
        seen.add(r[0])
        out.append(dict(zip(CORPUS_COLS, r)))
    return out


def read_summary(inp, want_taxids=None, want_accs=None):
    """GenBank assembly_summary.txt. The '#assembly_accession' header names the
    columns; every data line must have exactly that many fields (a short or
    long line is a truncated or wrong file, never skipped). Returns rows for
    the wanted species_taxids and/or accessions."""
    path = inp.path
    hdr, out, n_data = None, [], 0
    for n, line in enumerate(inp.lines(), 1):
        line = line.rstrip("\n")
        if line.startswith("#"):
            if line.startswith("#assembly_accession"):
                if hdr is not None:
                    raise Fail("%s line %d: second header line" % (path, n))
                hdr = line[1:].split("\t")
                miss = [c for c in SUMMARY_NEED if c not in hdr]
                if miss:
                    raise Fail("%s: header lacks %s" % (path, miss))
                ix = {c: hdr.index(c) for c in SUMMARY_NEED}
            continue
        if hdr is None:
            raise Fail("%s line %d: data before the #assembly_accession header" % (path, n))
        f = line.split("\t")
        if len(f) != len(hdr):
            raise Fail("%s line %d: %d fields, header has %d" % (path, n, len(f), len(hdr)))
        n_data += 1
        acc = f[ix["assembly_accession"]]
        if ((want_taxids is not None and f[ix["species_taxid"]] in want_taxids)
                or (want_accs is not None and acc in want_accs)):
            out.append({c: f[ix[c]] for c in SUMMARY_NEED})
    if hdr is None or n_data == 0:
        raise Fail("%s has no header or no data lines" % path)
    return out, n_data


def list_dir(d):
    """Names in d via scandir; a transient Lustre error is retried, a
    persistent one is fatal (never a short listing)."""
    last = None
    for t in range(LIST_TRIES):
        try:
            with os.scandir(d) as it:
                return sorted((e.name, e.is_dir(follow_symlinks=False)) for e in it)
        except FileNotFoundError:
            raise Fail("source directory %s does not exist" % d)
        except OSError as e:
            last = e
            time.sleep(0.5 * (t + 1))
    raise Fail("cannot list %s after %d tries: %s" % (d, LIST_TRIES, last))


# --------------------------------------------------------------- run frame
class Run:
    """New output directory, staged outputs, record written last."""

    def __init__(self, args, cmd):
        if not re.fullmatch(r"[0-9a-f]{40}", os.environ.get("GIT_COMMIT", "")):
            raise Fail("GIT_COMMIT is not set to a 40-hex commit (G-4): run as "
                       "GIT_COMMIT=$(git rev-parse HEAD) ... or sbatch --export=ALL,GIT_COMMIT=...")
        self.cmd, self.args = cmd, args
        git_provenance()                    # fails now, before any work, if unmeasurable
        self.out = os.path.abspath(args.out)
        self.tmp = os.path.join(self.out, ".tmp")
        self.inputs, self.outputs = {}, []
        self.created = False                # True once this run made self.out

    def start(self):
        """Create the output directory (refused if it exists) and .tmp. Any
        failure after the directory exists is cleaned up by fail(), which
        main() calls; a directory this run did not create is never touched."""
        if os.path.lexists(self.out):
            raise Fail("output directory %s already exists; each run needs a new one" % self.out)
        try:
            os.makedirs(self.out)
        except FileExistsError:
            raise Fail("output directory %s already exists; each run needs a new one" % self.out)
        self.created = True
        os.makedirs(self.tmp)

    def input(self, key, inp):
        """Record an Input read to the end (its own bytes and sha256)."""
        if inp.rec is None:
            raise Fail("%s %s was not read to the end" % (inp.what, inp.path))
        self.inputs[key] = inp.rec

    def path(self, name):
        self.outputs.append(name)
        return os.path.join(self.tmp, name)

    def publish(self, counts, extra=None):
        """Digests of every output first, then the moves, then the record.
        A failure anywhere in here is cleaned up by fail()."""
        shas = {}
        for name in self.outputs:
            src = os.path.join(self.tmp, name)
            shas[name] = {"bytes": os.path.getsize(src), "sha256": sha256_file(src)}
        rec = {"status": "complete", "command": self.cmd, "argv": sys.argv,
               "params": {k: v for k, v in vars(self.args).items() if k != "func"},
               "inputs": self.inputs, "counts": counts, "outputs": shas,
               "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
               "host": socket.gethostname(), "python": sys.version.split()[0],
               "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
               "slurm_array_task_id": os.environ.get("SLURM_ARRAY_TASK_ID")}
        rec.update(extra or {})
        rec.update(git_provenance())
        tmp = os.path.join(self.out, self.cmd + ".json.tmp")
        with open(tmp, "w") as fh:
            json.dump(rec, fh, indent=1, sort_keys=True)
            fh.write("\n")
        for name in self.outputs:
            os.replace(os.path.join(self.tmp, name), os.path.join(self.out, name))
        os.rmdir(self.tmp)
        os.replace(tmp, os.path.join(self.out, self.cmd + ".json"))     # last

    def fail(self, reason):
        """No output survives a failure: staged files, outputs already moved,
        and the record (temporary or final) are all removed. Nothing is done
        to an output directory this run did not create."""
        if not self.created:
            return
        shutil.rmtree(self.tmp, ignore_errors=True)
        for name in self.outputs + [self.cmd + ".json.tmp", self.cmd + ".json"]:
            try:
                os.remove(os.path.join(self.out, name))
            except FileNotFoundError:
                pass
            except OSError as e:
                log("cannot remove %s after failure: %s" % (name, e))
        rec = {"status": "failed", "command": self.cmd, "reason": reason, "argv": sys.argv,
               "host": socket.gethostname(), "slurm_job_id": os.environ.get("SLURM_JOB_ID")}
        try:
            rec.update(git_provenance())
        except Fail as e:
            rec["git"] = str(e)
        try:
            with open(os.path.join(self.out, "FAILED.json"), "w") as fh:
                json.dump(rec, fh, indent=1, sort_keys=True)
                fh.write("\n")
        except OSError:
            pass


def write_tsv(path, cols, rows):
    with open(path, "w", newline="") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            vals = [str(r[c]) for c in cols]
            if any(("\t" in v or "\n" in v) for v in vals):
                raise Fail("field with a tab or newline in %s: %r" % (path, vals))
            fh.write("\t".join(vals) + "\n")


# ------------------------------------------------------------------- build
def cmd_build(args, run):
    inp = Input(args.sources, "sources")
    sources = read_sources(inp)
    run.input("sources", inp)
    rows, where, per = [], {}, collections.OrderedDict()
    n_links = 0
    for s in sources:
        d = s["source_dir"]
        names = list_dir(d)
        n = 0
        for name, is_dir in names:
            if not name.endswith(".fna"):
                continue                    # .fai, _plasmid_only/, ...: not genomes
            path = d + "/" + name
            if is_dir:
                raise Fail("%s is a directory, not a genome file" % path)
            acc = name[:-4]
            acc_parts(acc)
            try:
                st = os.stat(path)         # follows a symlink: a broken one fails here
            except FileNotFoundError:
                raise Fail("%s is a broken symlink" % path)
            if not os.path.isfile(path):
                raise Fail("%s is not a regular file" % path)
            if st.st_size == 0:
                raise Fail("%s is empty" % path)
            if acc in where:
                raise Fail("assembly %s is in two sources: %s and %s" % (acc, where[acc], d))
            where[acc] = d
            n_links += os.path.islink(path)
            rows.append({"assembly": acc, "species": s["species"], "fna_path": path})
            n += 1
        if n != s["expected_n"]:
            raise Fail("%s holds %d genome files, config declares %d"
                       % (d, n, s["expected_n"]))
        per[s["species"]] = n
        log("%-15s %6d  %s" % (s["species"], n, d))
    rows.sort(key=lambda r: r["assembly"])
    out = run.path("corpus.tsv")
    write_tsv(out, CORPUS_COLS, rows)
    got = sha256_file(out)
    if args.expect_sha256 is not None:
        if not HEX64.fullmatch(args.expect_sha256):
            raise Fail("--expect-sha256 must be 64 lowercase hex")
        if got != args.expect_sha256:
            raise Fail("rebuilt corpus has sha256 %s, expected %s" % (got, args.expect_sha256))
    log("corpus: %d assemblies from %d sources, sha256 %s" % (len(rows), len(sources), got))
    run.publish({"assemblies": len(rows), "sources": len(sources), "symlinks": n_links,
                 "per_species": per},
                {"corpus_sha256": got, "expect_sha256_matched": args.expect_sha256 is not None})


# ------------------------------------------------------------------- pairs
def read_genome(path):
    """One complete read of a FASTA: file bytes and sha256, and its sequence
    content (module docstring). The bytes read must equal the file's size at
    open and after the read: an early EOF (a short read on a shared file
    system) or a file changing under us is fatal, never a shorter genome."""
    fh_sha, digests, h, n_bases, n_read = hashlib.sha256(), [], None, 0, 0
    strip = bytes.maketrans(b"", b"")
    with open(path, "rb") as fh:
        size0 = os.fstat(fh.fileno()).st_size
        for line in fh:
            fh_sha.update(line)
            n_read += len(line)
            if line.startswith(b">"):
                if h is not None:
                    digests.append(h.hexdigest())
                h = hashlib.sha256()
                continue
            s = line.translate(strip, ASCII_WS)
            if not s:
                continue
            if h is None:
                raise Fail("%s: sequence before the first header" % path)
            s = s.upper()
            h.update(s)
            n_bases += len(s)
        size1 = os.fstat(fh.fileno()).st_size
    if n_read != size0 or size1 != size0:
        raise Fail("%s: read %d bytes, file size %d at open and %d after (short read or "
                   "file changed)" % (path, n_read, size0, size1))
    if h is None:
        raise Fail("%s has no FASTA record" % path)
    digests.append(h.hexdigest())
    return {"bytes": n_read, "file_sha256": fh_sha.hexdigest(),
            "seq_sha256": hashlib.sha256("\n".join(sorted(digests)).encode()).hexdigest(),
            "records": len(digests), "bases": n_bases}


def _read_worker(item):
    acc, path = item
    try:
        return acc, read_genome(path), None
    except (Fail, OSError) as e:
        return acc, None, "%s: %s" % (type(e).__name__, e)


def read_all(items, workers):
    out = {}
    with multiprocessing.get_context("fork").Pool(workers) as pool:
        for i, (acc, sig, err) in enumerate(pool.imap_unordered(_read_worker, items, 8), 1):
            if err:
                raise Fail("cannot read %s (%s)" % (acc, err))
            out[acc] = sig
            if i % 2000 == 0:
                log("  %d/%d genomes read" % (i, len(items)))
    if len(out) != len(items):
        raise Fail("read %d of %d genomes" % (len(out), len(items)))
    return out


def canonical_of(members):
    """GCA before GCF, then the highest version: the GenBank record of the
    newest version stands for the assembly."""
    return max(members, key=lambda a: (acc_parts(a)[0] == "GCA", acc_parts(a)[2], a))


def cmd_pairs(args, run):
    inp = Input(args.corpus, "corpus", args.corpus_sha256)
    corpus = read_corpus(inp)
    run.input("corpus", inp)
    path_of = {r["assembly"]: r["fna_path"] for r in corpus}
    sp_of = {r["assembly"]: r["species"] for r in corpus}

    by_id = collections.defaultdict(list)
    for a in path_of:
        by_id[acc_parts(a)[1]].append(a)
    groups = sorted(sorted(v) for v in by_id.values() if len(v) > 1)

    gca = {a for g in groups for a in g if acc_parts(a)[0] == "GCA"}
    inp = Input(args.summary, "summary", args.summary_sha256)
    srows, _ = read_summary(inp, want_accs=gca)
    run.input("summary", inp)
    ncbi = {r["assembly_accession"]: r for r in srows}

    log("reading all %d genomes (%d same-id groups)" % (len(corpus), len(groups)))
    sig = read_all([(r["assembly"], r["fna_path"]) for r in corpus], args.workers)
    content = [dict(assembly=r["assembly"], species=r["species"], fna_path=r["fna_path"],
                    realpath=os.path.realpath(r["fna_path"]), **sig[r["assembly"]])
               for r in corpus]
    write_tsv(run.path("corpus_content.tsv"),
              ["assembly", "species", "fna_path", "realpath", "bytes", "file_sha256",
               "seq_sha256", "records", "bases"], content)

    dup_rows, rel_n = [], collections.Counter()
    for g in groups:
        can = canonical_of(g)
        for a in g:
            if a == can:
                continue
            if sig[a]["seq_sha256"] == sig[can]["seq_sha256"]:
                rel = "identical_sequence"
            elif acc_parts(a)[2] != acc_parts(can)[2]:
                rel = "different_version"
            else:
                rel = "differs"
            rel_n[rel] += 1
            nr = ncbi.get(can, {})
            dup_rows.append({
                "duplicate": a, "canonical": can, "species": sp_of[a], "relation": rel,
                "collapse": "yes" if rel != "differs" else "no",
                "dup_seq_sha256": sig[a]["seq_sha256"],
                "canonical_seq_sha256": sig[can]["seq_sha256"],
                "dup_records": sig[a]["records"], "canonical_records": sig[can]["records"],
                "dup_bases": sig[a]["bases"], "canonical_bases": sig[can]["bases"],
                "ncbi_paired_asm": nr.get("gbrs_paired_asm", "not_in_summary"),
                "ncbi_paired_comp": nr.get("paired_asm_comp", "not_in_summary")})
    write_tsv(run.path("corpus_duplicates.tsv"),
              ["duplicate", "canonical", "species", "relation", "collapse",
               "dup_seq_sha256", "canonical_seq_sha256", "dup_records", "canonical_records",
               "dup_bases", "canonical_bases", "ncbi_paired_asm", "ncbi_paired_comp"], dup_rows)

    by_seq = collections.defaultdict(list)
    for a in path_of:
        by_seq[sig[a]["seq_sha256"]].append(a)
    ident_rows, n_groups = [], 0
    for s_, members in sorted(by_seq.items(), key=lambda kv: sorted(kv[1])[0]):
        ids = {acc_parts(a)[1] for a in members}
        if len(ids) < 2:
            continue                        # one assembly, or one id (a same-id group)
        n_groups += 1
        for a in sorted(members):
            ident_rows.append({"group_id": "I%05d" % n_groups, "assembly": a,
                               "species": sp_of[a], "seq_sha256": s_,
                               "records": sig[a]["records"], "bases": sig[a]["bases"],
                               "group_size": len(members), "distinct_ids": len(ids)})
    write_tsv(run.path("corpus_identical_other_id.tsv"),
              ["group_id", "assembly", "species", "seq_sha256", "records", "bases",
               "group_size", "distinct_ids"], ident_rows)

    n_collapse = sum(r["collapse"] == "yes" for r in dup_rows)
    counts = {"corpus_assemblies": len(corpus), "genomes_read": len(sig),
              "bytes_read": sum(v["bytes"] for v in sig.values()),
              "same_id_groups": len(groups), "relation": dict(rel_n),
              "collapse_yes": n_collapse, "assemblies_after_collapse": len(corpus) - n_collapse,
              "identical_other_id_groups": n_groups,
              "identical_other_id_assemblies": len(ident_rows)}
    log(json.dumps(counts, sort_keys=True))
    run.publish(counts, {"corpus_sha256": args.corpus_sha256})


# -------------------------------------------------------------------- gaps
def cmd_gaps(args, run):
    if not DATE.fullmatch(args.cutoff):
        raise Fail("--cutoff must be YYYY-MM-DD, got %r" % args.cutoff)
    inp = Input(args.sources, "sources")
    sources = read_sources(inp)
    run.input("sources", inp)
    inp = Input(args.corpus, "corpus", args.corpus_sha256)
    corpus = read_corpus(inp)
    run.input("corpus", inp)
    species = {s["species"]: s for s in sources}
    unknown = sorted({r["species"] for r in corpus} - set(species))
    if unknown:
        raise Fail("corpus species %s are not declared in %s" % (unknown, args.sources))
    tax_sp = {s["species_taxid"]: s["species"] for s in sources}
    if len(tax_sp) != len(sources):
        raise Fail("two sources share a species_taxid")

    have = collections.defaultdict(dict)        # species -> id -> corpus accession(s)
    n_corpus, n_gcf = collections.Counter(), collections.Counter()
    for r in corpus:
        p, i, _ = acc_parts(r["assembly"])
        have[r["species"]].setdefault(i, []).append(r["assembly"])
        n_corpus[r["species"]] += 1
        n_gcf[r["species"]] += p == "GCF"
    plas = collections.defaultdict(dict)        # species -> id -> [accessions], every version
    plas_rows = []
    for s in sources:
        pd = os.path.join(s["source_dir"], PLASMID_DIR)
        # only a path with no directory entry at all means "no set-aside
        # genomes"; an entry that exists but cannot be followed or checked
        # (a broken symlink, no permission) is fatal
        try:
            os.lstat(pd)
        except FileNotFoundError:
            continue                        # no set-aside genomes for this species
        except OSError as e:
            raise Fail("cannot check %s: %s" % (pd, e))
        try:
            st = os.stat(pd)
        except OSError as e:
            raise Fail("%s exists but cannot be checked: %s" % (pd, e))
        if not stat.S_ISDIR(st.st_mode):
            raise Fail("%s exists but is not a directory" % pd)
        for name, is_dir in list_dir(pd):
            if not name.endswith(".fna"):
                continue
            path = pd + "/" + name
            if is_dir:
                raise Fail("%s is a directory, not a genome file" % path)
            acc = name[:-4]
            acc_parts(acc)
            try:
                size = os.stat(path).st_size
            except OSError as e:
                raise Fail("cannot stat %s: %s" % (path, e))
            plas[s["species"]].setdefault(acc_parts(acc)[1], []).append(acc)
            plas_rows.append({"species": s["species"], "accession": acc, "path": path,
                              "bytes": size})
    # the set-aside listing is an input of this report: published with it
    write_tsv(run.path("corpus_plasmid_only.tsv"), ["species", "accession", "path", "bytes"],
              sorted(plas_rows, key=lambda r: (r["species"], r["accession"])))

    inp = Input(args.summary, "summary", args.summary_sha256)
    srows, n_data = read_summary(inp, want_taxids=set(tax_sp))
    run.input("summary", inp)
    usable = collections.defaultdict(dict)       # species -> id -> row
    for r in srows:
        if r["assembly_level"] not in USABLE_LEVELS:
            continue
        if not DATE.fullmatch(r["seq_rel_date"]):
            raise Fail("summary: %s has seq_rel_date %r (expected YYYY-MM-DD)"
                       % (r["assembly_accession"], r["seq_rel_date"]))
        p, i, _ = acc_parts(r["assembly_accession"])
        if p != "GCA":
            raise Fail("summary row %s is not a GenBank (GCA) accession; use the "
                       "GenBank summary" % r["assembly_accession"])
        sp = tax_sp[r["species_taxid"]]
        if i in usable[sp]:
            raise Fail("summary lists id %s twice for %s" % (i, sp))
        usable[sp][i] = r

    per, gap = [], []
    for s in sources:
        sp = s["species"]
        u, h, pl = usable[sp], have[sp], plas[sp]
        covered = set(h) | set(pl)
        n_pl = sum(len(v) for v in pl.values())     # listed files, every version counted
        c = collections.Counter()
        for i, r in sorted(u.items()):
            at = r["seq_rel_date"] <= args.cutoff
            c["usable_at_cutoff"] += at
            if i in covered:
                # the GenBank accessions held for this id (corpus and set-aside,
                # every version); a mismatch only if none is the summary's
                held = sorted([a for a in h.get(i, []) if a.startswith("GCA")]
                              + [a for a in pl.get(i, []) if a.startswith("GCA")],
                              key=acc_parts)                # by version number
                if held and r["assembly_accession"] not in held:
                    c["version_mismatch"] += 1
                    gap.append(dict(r, species=sp, status="version_mismatch",
                                    corpus_accession=",".join(held)))
                continue
            st = "missing_at_cutoff" if at else "released_after_cutoff"
            c[st] += 1
            gap.append(dict(r, species=sp, status=st, corpus_accession="."))
        for i in sorted(covered - set(u)):
            acc = (h.get(i) or sorted(pl[i], key=acc_parts))[0]
            c["not_in_usable_set"] += 1
            gap.append({"species": sp, "assembly_accession": acc, "status": "not_in_usable_set",
                        "seq_rel_date": ".", "assembly_level": ".", "ftp_path": ".",
                        "corpus_accession": acc})
        per.append({"species": sp, "species_taxid": s["species_taxid"],
                    "corpus_assemblies": n_corpus[sp], "corpus_gcf": n_gcf[sp],
                    "plasmid_only_set_aside": n_pl, "usable_total": len(u),
                    "usable_at_cutoff": c["usable_at_cutoff"],
                    "missing_at_cutoff": c["missing_at_cutoff"],
                    "released_after_cutoff": c["released_after_cutoff"],
                    "not_in_usable_set": c["not_in_usable_set"],
                    "version_mismatch": c["version_mismatch"]})
    cols = ["species", "species_taxid", "corpus_assemblies", "corpus_gcf",
            "plasmid_only_set_aside", "usable_total", "usable_at_cutoff",
            "missing_at_cutoff", "released_after_cutoff", "not_in_usable_set",
            "version_mismatch"]
    write_tsv(run.path("corpus_gaps.tsv"), cols, per)
    gap.sort(key=lambda r: (r["species"], r["status"], r["assembly_accession"]))
    write_tsv(run.path("corpus_gap_list.tsv"),
              ["species", "accession", "status", "seq_rel_date", "assembly_level",
               "ftp_path", "corpus_accession"],
              [dict(r, accession=r["assembly_accession"]) for r in gap])
    # the download queue: missing_at_cutoff rows that have an NCBI ftp_path
    # (a row without one stays in the gap list and is counted, never queued)
    queue, no_ftp = [], 0
    for r in gap:
        if r["status"] != "missing_at_cutoff":
            continue
        base = os.path.basename(r["ftp_path"].rstrip("/"))
        acc = r["assembly_accession"]
        if not FTP.fullmatch(r["ftp_path"]) or not (base == acc or base.startswith(acc + "_")):
            no_ftp += 1
            continue
        queue.append(dict(r, accession=r["assembly_accession"]))
    write_tsv(run.path("corpus_download_queue.tsv"),
              ["species", "accession", "seq_rel_date", "assembly_level", "ftp_path"], queue)
    tot = {k: sum(r[k] for r in per) for k in cols[2:]}
    log(json.dumps(tot, sort_keys=True))
    run.publish({"summary_data_lines": n_data, "totals": tot, "queued": len(queue),
                 "missing_without_ftp_path": no_ftp, "plasmid_only_listed": len(plas_rows)},
                {"cutoff": args.cutoff})


# -------------------------------------------------------------------- main
def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build")
    b.add_argument("--sources", required=True)
    b.add_argument("--expect-sha256", default=None)
    b.add_argument("--out", required=True)
    p = sub.add_parser("pairs")
    for k in ("corpus", "summary"):
        p.add_argument("--" + k, required=True)
        p.add_argument("--%s-sha256" % k, required=True)
    p.add_argument("--workers", type=int, default=8)
    p.add_argument("--out", required=True)
    g = sub.add_parser("gaps")
    g.add_argument("--sources", required=True)
    for k in ("corpus", "summary"):
        g.add_argument("--" + k, required=True)
        g.add_argument("--%s-sha256" % k, required=True)
    g.add_argument("--cutoff", required=True, help="YYYY-MM-DD; seq_rel_date <= cutoff")
    g.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    if a.cmd == "pairs" and a.workers < 1:
        ap.error("--workers must be >= 1")
    return a


def main(argv=None):
    args = parse_args(argv)
    run = None
    try:
        run = Run(args, args.cmd)
        run.start()
        {"build": cmd_build, "pairs": cmd_pairs, "gaps": cmd_gaps}[args.cmd](args, run)
        return 0
    except Fail as e:
        log("FATAL: %s" % e)
        if run:
            run.fail(str(e))
        return 2
    except Exception as e:          # OSError, bugs: same contract, never a partial result
        import traceback
        traceback.print_exc()
        log("FATAL: %s: %s" % (type(e).__name__, e))
        if run:
            run.fail("%s: %s" % (type(e).__name__, e))
        return 2


if __name__ == "__main__":
    sys.exit(main())
