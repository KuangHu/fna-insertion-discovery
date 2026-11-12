#!/usr/bin/env python3
"""Tests for scripts/248_trusted_manifest.py (SPEC item 0, V0-6).

Small synthetic fixtures in a temporary directory; no real data, no BLAST
(the blast bin is a stub that only answers -version).  Run on the login node:
    python3 -I tests/test_trusted_manifest.py
"""
import argparse
import errno
import functools
import hashlib
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(HERE, "..", "scripts", "248_trusted_manifest.py")


def load_248():
    spec = importlib.util.spec_from_file_location("m248", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    sys.modules["m248"] = m        # Pool pickles the worker function by module name
    return m


def sha(b):
    return hashlib.sha256(b).hexdigest()


class Fixture:
    """5 genomes in 2 shards (sizes 3 + 2), 3 queries."""

    def __init__(self, root):
        self.root = root
        self.g = os.path.join(root, "genomes")
        os.makedirs(self.g)
        self.genomes = {
            "GCA_1.1": b">c1 desc\nACGTACGT\nAC\n>c2\nGGGG\n",
            "GCA_2.1": b">c1\nAAAA\n",
            "GCA_3.1": b">x\tTAB\nCCCCCC\n",
            "GCA_4.1": b">p1\nAC\n\n>p2\nGT\n",
            "GCA_5.1": b">only\nACGTN\n",
        }
        self.queries_fna = b">Q_aa\nACGTACGTAC\n>Q_bb\nGGGGCCCC\n>Q_cc extra words\nAT\nAT\n"
        self.queries_tsv = (b"query_id\tbag_id\tseq_key\tlen_nt\n"
                            b"Q_aa\tB1\taa\t10\nQ_bb\tB1\tbb\t8\nQ_cc\tB2\tcc\t4\n")
        self.order = list(self.genomes)
        self.blast = os.path.join(root, "bl")
        os.makedirs(self.blast)
        self.write_blast("2.17.0")
        # a throwaway git repo standing in for the project repo (G-4 / code check)
        self.repo = os.path.join(root, "repo")
        os.makedirs(os.path.join(self.repo, "code"))
        for name in ("a.sh", "b.py"):
            with open(os.path.join(self.repo, "code", name), "w") as f:
                f.write(f"# {name}\n")
        git = lambda *a: subprocess.run(["git", "-C", self.repo, *a], check=True, capture_output=True)
        git("init", "-q"); git("config", "user.email", "t@t"); git("config", "user.name", "t")
        git("add", "-A"); git("commit", "-qm", "c")
        self.commit = subprocess.run(["git", "-C", self.repo, "rev-parse", "HEAD"],
                                     capture_output=True, text=True, check=True).stdout.strip()
        self.code_files = ["code/a.sh", "code/b.py"]
        # tool self-test outputs as V0_prep.sh leaves them
        self.tt = os.path.join(root, "tooltests")
        os.makedirs(self.tt)
        H = "a" * 64
        for name, text in (("outfmt7_test.verdict", "V0 TEST1 PASS: all 3\n"),
                           ("v0b_dbtest.log",
                            f"git_commit_submit={'1'*40} git_commit_runtime={'1'*40} git_dirty=false\n"
                            f"script sha256: {H}\ninput: 3 contigs; x.fna sha256 {H}\nV0B PASS: ok\n"),
                           ("commands.txt",
                            f"# V0 job 1; test genome /g.fna (sha256 {H})\n"
                            f"# TEST 1 test inputs: t_q.fna sha256 {H}, t_db.fna sha256 {H}\n"
                            "blastn -task dc-megablast -evalue 1e-30\n")):
            with open(os.path.join(self.tt, name), "w") as f:
                f.write(text)
        self.write()

    def write_blast(self, version):
        for t in ("blastn", "makeblastdb", "blastdbcmd"):
            p = os.path.join(self.blast, t)
            with open(p, "w") as f:
                f.write(f"#!/bin/sh\necho '{t}: {version}+'\necho 'Package: blast {version}, build x'\n")
            os.chmod(p, 0o755)

    def write(self, corpus_rows=None):
        for acc, data in self.genomes.items():
            with open(os.path.join(self.g, acc + ".fna"), "wb") as f:
                f.write(data)
        with open(os.path.join(self.root, "q.fna"), "wb") as f:
            f.write(self.queries_fna)
        with open(os.path.join(self.root, "q.tsv"), "wb") as f:
            f.write(self.queries_tsv)
        rows = corpus_rows if corpus_rows is not None else self.order
        with open(os.path.join(self.root, "corpus.tsv"), "w") as f:
            f.write("assembly\tspecies\tfna_path\n")
            for acc in rows:
                f.write(f"{acc}\tsp\t{os.path.join(self.g, acc + '.fna')}\n")
        ids = sorted(l[1:].split()[0] for l in self.queries_fna.split(b"\n") if l.startswith(b">"))
        self.expect = {
            "all_queries_fna_sha256": sha(self.queries_fna),
            "query_ids_sorted_sha256": sha(b"".join(i + b"\n" for i in ids)),
            "n_queries": len(ids),
            "all_queries_fna_bytes": len(self.queries_fna),
            "n_genomes": len(rows),
            "n_shards": 2,
            "shard_sizes": "1x3,1x2",
            "blast_version": "2.17.0",
        }
        self.write_expect()

    def write_expect(self):
        with open(os.path.join(self.root, "expect.json"), "w") as f:
            json.dump(self.expect, f)

    def build(self, out="run", env_commit="HEAD"):
        env = dict(os.environ)
        env.pop("GIT_COMMIT", None)
        if env_commit == "HEAD":
            env_commit = self.commit
        if env_commit is not None:
            env["GIT_COMMIT"] = env_commit
        cmd = [sys.executable, "-I", SCRIPT, "build",
               "--queries-fna", os.path.join(self.root, "q.fna"),
               "--queries-tsv", os.path.join(self.root, "q.tsv"),
               "--corpus", os.path.join(self.root, "corpus.tsv"),
               "--expect", os.path.join(self.root, "expect.json"),
               "--blast-bin", self.blast, "--repo", self.repo, "--threads", "2",
               "--tooltests", self.tt, "--code-files", *self.code_files,
               "--out", os.path.join(self.root, out)]
        return subprocess.run(cmd, capture_output=True, text=True, env=env)

    def check(self, out="run"):
        cmd = [sys.executable, "-I", SCRIPT, "check", "--out", os.path.join(self.root, out),
               "--corpus", os.path.join(self.root, "corpus.tsv"),
               "--expect", os.path.join(self.root, "expect.json")]
        return subprocess.run(cmd, capture_output=True, text=True)


class T(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="t248.", dir=os.environ.get("TMPDIR"))
        self.f = Fixture(self.tmp)

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def out(self, name="run"):
        return os.path.join(self.tmp, name)

    def assertFails(self, r, frag):
        self.assertNotEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertIn(frag, r.stderr, r.stderr)

    def assertBuildFails(self, frag):
        r = self.f.build()
        self.assertFails(r, frag)
        self.assertFalse(os.path.exists(os.path.join(self.out(), "v0_result.json")))
        self.assertTrue(os.path.exists(os.path.join(self.out(), "FAILED")))

    def read(self, name):
        with open(os.path.join(self.out(), name)) as f:
            return f.read()

    # ---------------- positive
    def test_build_ok(self):
        r = self.f.build()
        self.assertEqual(r.returncode, 0, r.stderr)
        res = json.loads(self.read("v0_result.json"))
        self.assertTrue(res["verified"])
        self.assertEqual(res["counts"], {"queries": 3, "genomes": 5, "contigs": 7, "bases": 33})
        self.assertEqual(res["git_commit_submit"], self.f.commit)
        self.assertIs(res["git_dirty"], False)
        self.assertEqual(set(res["code_files_sha256_at_commit"]), set(self.f.code_files))
        self.assertIn("dc-megablast", res["tooltests"]["commands.txt"]["text"])
        self.assertEqual(set(res["tooltests"]), {"outfmt7_test.verdict", "v0b_dbtest.log", "commands.txt"})
        for k in ("git_commit_runtime", "git_dirty", "host", "tools", "argv", "inputs"):
            self.assertIn(k, res)
        contigs = self.read("contigs.tsv").splitlines()
        self.assertIn("0\tGCA_1.1\tGCA_1.1|c1\t10", contigs)       # header token only; 2 lines summed
        self.assertIn("0\tGCA_3.1\tGCA_3.1|x\t6", contigs)         # tab ends the id (awk FS)
        self.assertIn("1\tGCA_4.1\tGCA_4.1|p1\t2", contigs)        # blank line contributes 0
        genomes = self.read("genomes.tsv").splitlines()
        self.assertEqual([l.split("\t")[0] for l in genomes[1:]], ["0", "1", "0", "1", "0"])
        g1 = genomes[1].split("\t")
        self.assertEqual(g1[4], sha(self.f.genomes["GCA_1.1"]))
        self.assertEqual(g1[7], sha(b"c1\nc2\n"))
        t = dict(l.split("\t")[:2] for l in self.read("isref_all_trusted.tsv").splitlines()[1:])
        self.assertEqual(t["genome_manifest_sha256"], sha(open(os.path.join(self.out(), "genome_manifest.index"), "rb").read()))
        self.assertEqual(t["all_queries_tsv_sha256"], sha(self.f.queries_tsv))
        self.assertEqual(self.f.check().returncode, 0)

    def test_never_overwrites(self):
        os.makedirs(self.out())
        r = self.f.build()
        self.assertFails(r, "out_exists")
        self.assertEqual(os.listdir(self.out()), [])

    def test_git_commit_required(self):
        r = self.f.build(env_commit=None)
        self.assertFails(r, "git_commit_submit_unset")

    # ---------------- G-4 / V0-5: provenance of the code that ran
    def test_untracked_file_makes_dirty(self):
        open(os.path.join(self.f.repo, "code", "new_untracked.py"), "w").write("x\n")
        r = self.f.build()
        self.assertEqual(r.returncode, 0, r.stderr)
        res = json.loads(self.read("v0_result.json"))
        self.assertIs(res["git_dirty"], True)
        self.assertIn("?? code/new_untracked.py", res["git_status_porcelain"])

    def test_code_modified_refused(self):
        open(os.path.join(self.f.repo, "code", "b.py"), "a").write("changed\n")
        self.assertBuildFails("code_differs_from_commit: code/b.py")

    def test_code_untracked_refused(self):
        open(os.path.join(self.f.repo, "code", "c.sh"), "w").write("x\n")
        self.f.code_files.append("code/c.sh")
        self.assertBuildFails("code_not_committed: code/c.sh")

    def test_blast_version_pinned(self):
        self.f.write_blast("2.16.0")
        self.assertBuildFails("blast_version_mismatch: blastn 2.16.0 != 2.17.0")

    def test_tooltest_digest_missing(self):
        p = os.path.join(self.f.tt, "commands.txt")
        data = open(p).read().replace("t_db.fna sha256 " + "a" * 64, "t_db.fna sha256 ")
        open(p, "w").write(data)
        self.assertBuildFails("tooltest_digest_missing: commands.txt")

    def test_tooltest_missing(self):
        os.remove(os.path.join(self.f.tt, "commands.txt"))
        self.assertBuildFails("tooltest_missing: commands.txt")

    def test_tooltest_failed(self):
        open(os.path.join(self.f.tt, "outfmt7_test.verdict"), "w").write("V0 TEST1 FAIL\n")
        self.assertBuildFails("tooltest_not_passed: outfmt7_test.verdict")

    # ---------------- V0-6 negatives: missing file, duplicate id, intersection, count, length sum
    def test_missing_genome_file(self):
        os.remove(os.path.join(self.f.g, "GCA_4.1.fna"))
        self.assertBuildFails("genome_missing: GCA_4.1")

    def test_duplicate_query_id(self):
        self.f.queries_fna = b">Q_aa\nACGTACGTAC\n>Q_aa\nGGGGCCCC\n>Q_cc\nATAT\n"
        self.f.write()
        self.assertBuildFails("duplicate_query_id")

    def test_duplicate_contig_id(self):
        self.f.genomes["GCA_2.1"] = b">c1\nAAAA\n>c1 again\nCC\n"
        self.f.write()
        self.assertBuildFails("duplicate_contig_id: GCA_2.1")

    def test_shards_intersect_duplicate_assembly(self):
        self.f.write(corpus_rows=["GCA_1.1", "GCA_2.1", "GCA_3.1", "GCA_4.1", "GCA_1.1"])
        self.assertBuildFails("duplicate_assembly_in_corpus")

    def test_query_count_literal_mismatch(self):
        self.f.expect["n_queries"] = 4
        self.f.write_expect()
        self.assertBuildFails("literal_mismatch: n_queries")

    def test_fna_sha_literal_mismatch(self):
        self.f.expect["all_queries_fna_sha256"] = "f" * 64
        self.f.write_expect()
        self.assertBuildFails("literal_mismatch: all_queries_fna_sha256")

    def test_genome_count_mismatch(self):
        self.f.expect.update(n_genomes=6, shard_sizes="1x3,1x3")
        self.f.write_expect()
        self.assertBuildFails("genome_count")

    # ---------------- G-7 types of the expected values
    def _expect_type_fails(self, key, value, frag="expect_type"):
        self.f.expect[key] = value
        self.f.write_expect()
        self.assertBuildFails(frag)
        self.assertEqual(self.f.check().returncode != 0, True)

    def test_expect_float_count(self):
        self._expect_type_fails("n_queries", 3.0, "expect_type: n_queries")

    def test_expect_bool_count(self):
        self._expect_type_fails("n_shards", True, "expect_type: n_shards")

    def test_expect_string_count(self):
        self._expect_type_fails("all_queries_fna_bytes", "1", "expect_type: all_queries_fna_bytes")

    def test_expect_uppercase_sha(self):
        self._expect_type_fails("all_queries_fna_sha256", self.f.expect["all_queries_fna_sha256"].upper(),
                                "expect_type: all_queries_fna_sha256")

    def test_expect_unknown_key(self):
        self._expect_type_fails("n_quereis", 3, "expect_file_unknown_keys")

    def test_expect_sizes_inconsistent(self):
        self._expect_type_fails("shard_sizes", "1x3,1x3", "expect_inconsistent")

    def test_shard_size_mismatch(self):
        self.f.expect["shard_sizes"] = "1x4,1x1"   # adds up to 2 shards / 5 genomes, wrong split
        self.f.write_expect()
        self.assertBuildFails("shard_sizes")

    def test_contig_length_sum_vs_bases(self):
        self.assertEqual(self.f.build().returncode, 0)
        # tamper one contig length and re-seal the index, so only the sum check can catch it
        p = os.path.join(self.out(), "contigs.tsv")
        data = open(p).read().replace("GCA_1.1|c2\t4", "GCA_1.1|c2\t5")
        open(p, "w").write(data)
        g = open(os.path.join(self.out(), "genomes.tsv"), "rb").read()
        open(os.path.join(self.out(), "genome_manifest.index"), "w").write(
            f"{sha(g)}  genomes.tsv\n{sha(data.encode())}  contigs.tsv\n")
        self.assertFails(self.f.check(), "contig_length_sum: GCA_1.1")

    def test_tamper_without_reseal_caught_by_index(self):
        self.assertEqual(self.f.build().returncode, 0)
        p = os.path.join(self.out(), "contigs.tsv")
        data = open(p).read().replace("GCA_1.1|c2\t4", "GCA_1.1|c2\t5")
        open(p, "w").write(data)
        self.assertFails(self.f.check(), "index_sha_mismatch: contigs.tsv")

    def test_contig_count_vs_genomes(self):
        self.assertEqual(self.f.build().returncode, 0)
        p = os.path.join(self.out(), "contigs.tsv")
        lines = open(p).read().splitlines(True)
        data = "".join(l for l in lines if "GCA_1.1|c2" not in l)
        open(p, "w").write(data)
        g = open(os.path.join(self.out(), "genomes.tsv"), "rb").read()
        open(os.path.join(self.out(), "genome_manifest.index"), "w").write(
            f"{sha(g)}  genomes.tsv\n{sha(data.encode())}  contigs.tsv\n")
        self.assertFails(self.f.check(), "contig_r")   # count or order check fires

    def test_trusted_literal_tamper(self):
        self.assertEqual(self.f.build().returncode, 0)
        p = os.path.join(self.out(), "isref_all_trusted.tsv")
        data = open(p).read().replace("n_queries\t3\t", "n_queries\t03\t")
        open(p, "w").write(data)
        self.assertFails(self.f.check(), "type: n_queries")

    # ---------------- query tsv rules (V2-3(b))
    def test_tsv_len_mismatch(self):
        self.f.queries_tsv = self.f.queries_tsv.replace(b"bb\t8", b"bb\t9")
        self.f.write()
        self.assertBuildFails("queries_tsv_len_mismatch: Q_bb")

    def test_tsv_id_set(self):
        self.f.queries_tsv = b"query_id\tbag_id\tseq_key\tlen_nt\nQ_aa\tB1\taa\t10\nQ_bb\tB1\tbb\t8\n"
        self.f.write()
        self.assertBuildFails("queries_tsv_id_set")

    def test_tsv_empty_bag(self):
        self.f.queries_tsv = self.f.queries_tsv.replace(b"Q_cc\tB2", b"Q_cc\t")
        self.f.write()
        self.assertBuildFails("queries_tsv_empty_bag")

    # ---------------- strict FASTA rules
    def test_carriage_return(self):
        self.f.genomes["GCA_5.1"] = b">only\r\nACGT\r\n"
        self.f.write()
        self.assertBuildFails("carriage_return: GCA_5.1")

    def test_bad_sequence_char(self):
        self.f.genomes["GCA_5.1"] = b">only\nAC GT\n"
        self.f.write()
        self.assertBuildFails("bad_sequence_char: GCA_5.1")

    def test_zero_length_contig(self):
        self.f.genomes["GCA_5.1"] = b">a\n>b\nACGT\n"
        self.f.write()
        self.assertBuildFails("zero_length_contig: GCA_5.1")

    def test_empty_header(self):
        self.f.genomes["GCA_5.1"] = b">\nACGT\n"
        self.f.write()
        self.assertBuildFails("empty_id")

    def test_sequence_before_header(self):
        self.f.genomes["GCA_5.1"] = b"ACGT\n>a\nAC\n"
        self.f.write()
        self.assertBuildFails("sequence_before_header")

    def test_no_trailing_newline_ok(self):
        self.f.genomes["GCA_5.1"] = b">only\nACGTN"
        self.f.write()
        r = self.f.build()
        self.assertEqual(r.returncode, 0, r.stderr)


# ---------------------------------------------------------------- V0-7 I/O re-reads
class FaultyFile:
    """Wraps a real file; read() raises OSError(code) after the first block."""

    def __init__(self, f, code):
        self.f, self.code, self.n = f, code, 0

    def fileno(self):
        return self.f.fileno()

    def read(self, size=-1):
        self.n += 1
        if self.n > 1:
            raise OSError(self.code, os.strerror(self.code))
        return self.f.read(4)          # a short first block, then the error

    def __enter__(self):
        return self

    def __exit__(self, *a):
        self.f.close()


class FaultPlan:
    """Shared between forked workers through files: attempt k of path p fails
    with errno plan[p][k] (None = succeeds)."""

    def __init__(self, root, plan, where="read"):
        self.root, self.plan, self.where = root, plan, where
        os.makedirs(root, exist_ok=True)

    def attempt(self, path):
        cf = os.path.join(self.root, os.path.basename(path) + ".n")
        with open(cf, "a") as f:
            f.write("x")
        return os.path.getsize(cf) - 1

    def attempts(self, path):
        cf = os.path.join(self.root, os.path.basename(path) + ".n")
        return os.path.getsize(cf) if os.path.exists(cf) else 0

    def opener_for(self, path):
        k = self.attempt(path)
        codes = self.plan.get(os.path.basename(path), [])
        code = codes[k] if k < len(codes) else None

        def op(p, mode):
            if code is None:
                return open(p, mode)
            if self.where == "open":
                raise OSError(code, os.strerror(code))
            return FaultyFile(open(p, mode), code)
        return op


class TIO(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="t248io.", dir=os.environ.get("TMPDIR"))
        self.m = load_248()
        self.m.IO_BACKOFF_SECONDS = (0.0, 0.0, 0.0)
        self.orig_scan = self.m.scan_fasta
        self.p = os.path.join(self.tmp, "g.fna")
        with open(self.p, "wb") as f:
            f.write(b">c1 d\nACGTACGTAC\n>c2\nGG\n")
        self.sleeps = []
        self.counter = argparse.Namespace(value=0, get_lock=lambda: _NullLock())

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def plan(self, codes, where="read"):
        fp = FaultPlan(os.path.join(self.tmp, "plan"), {"g.fna": codes}, where)
        orig = self.orig_scan
        self.m.scan_fasta = lambda path, label: orig(path, label, opener=fp.opener_for(path))
        return fp

    def read(self):
        self.m.IO_BACKOFF_SECONDS = (10, 30)
        return self.m.read_genome(self.p, "G", self.counter, sleep=self.sleeps.append)

    def test_transient_error_recovers_with_confirming_read(self):
        fp = self.plan([errno.EIO])
        nbytes, sha256, recs, io = self.read()
        self.assertEqual(nbytes, os.path.getsize(self.p))
        self.assertEqual(sha256, sha(open(self.p, "rb").read()))
        self.assertEqual(recs, [(b"c1", 10), (b"c2", 2)])
        self.assertEqual(io, {"rereads": 2, "io_errors": ["EIO"]})   # failed read, good read, confirming read
        self.assertEqual(fp.attempts(self.p), 3)
        self.assertEqual(self.sleeps, [10])
        self.assertEqual(self.counter.value, 1)

    def test_error_logged_immediately(self):
        import io as _io
        self.plan([errno.EIO])
        ev, old = [], sys.stderr
        sys.stderr = buf = _io.StringIO()
        try:
            self.m.IO_BACKOFF_SECONDS = (10, 30)
            self.m.read_genome(self.p, "G", self.counter, sleep=self.sleeps.append, events=ev)
        finally:
            sys.stderr = old
        self.assertEqual([(d["assembly"], d["errno"], d["attempt"], d["node_total"]) for d in ev],
                         [("G", "EIO", 1, 1)])
        self.assertIn('IO_ERROR {', buf.getvalue())
        self.assertIn('"errno": "EIO"', buf.getvalue())

    def test_no_error_reads_once(self):
        fp = self.plan([])
        *_, io = self.read()
        self.assertEqual(io, {"rereads": 0, "io_errors": []})
        self.assertEqual(fp.attempts(self.p), 1)
        self.assertEqual(self.sleeps, [])

    def test_enotconn_at_open_and_eshutdown(self):
        fp = self.plan([errno.ENOTCONN, errno.ESHUTDOWN], where="open")
        *_, io = self.read()
        self.assertEqual(io, {"rereads": 3, "io_errors": ["ENOTCONN", "ESHUTDOWN"]})
        self.assertEqual(fp.attempts(self.p), 4)           # 2 failed, 1 complete, 1 confirming
        self.assertEqual(self.sleeps, [10, 30])            # increasing back-off

    def test_persistent_error_fails(self):
        fp = self.plan([errno.EIO] * 10)
        with self.assertRaises(self.m.GenomeReadFail) as cm:
            self.read()
        self.assertIn("genome_io_retries_exhausted: G", str(cm.exception))
        # error 3 leaves 1 read in the 4-read budget: no room for a complete
        # read plus a confirming read, so it fails without waiting again
        self.assertEqual(fp.attempts(self.p), 3)
        self.assertEqual(self.sleeps, [10, 30])
        self.assertEqual(cm.exception.io, {"rereads": 2, "io_errors": ["EIO"] * 3})

    def test_error_on_confirming_read_fails(self):
        # E, complete, E on the confirming read: 2 errors in, budget leaves 1 read
        fp = self.plan([errno.EIO, None, errno.EIO])
        with self.assertRaises(self.m.GenomeReadFail) as cm:
            self.read()
        self.assertIn("genome_io_retries_exhausted: G", str(cm.exception))
        self.assertEqual(fp.attempts(self.p), 3)

    def test_other_errno_not_retried(self):
        fp = self.plan([errno.EACCES])
        with self.assertRaises(self.m.GenomeReadFail) as cm:
            self.read()
        self.assertIn("genome_unreadable: G", str(cm.exception))
        self.assertEqual(fp.attempts(self.p), 1)
        self.assertEqual((self.sleeps, self.counter.value), ([], 0))

    def test_confirming_read_must_agree(self):
        calls = []

        def scan(path, label):
            calls.append(1)
            if len(calls) == 1:
                raise OSError(errno.EIO, "x")
            return (5, "a" * 64, [(b"c", 5)]) if len(calls) == 2 else (5, "b" * 64, [(b"c", 5)])
        self.m.scan_fasta = scan
        with self.assertRaises(self.m.GenomeReadFail) as cm:
            self.read()
        self.assertIn("genome_reread_mismatch: G", str(cm.exception))

    def test_node_limit(self):
        self.counter.value = 4
        self.plan([errno.EIO])
        with self.assertRaises(self.m.GenomeReadFail) as cm:
            self.read()
        self.assertIn("node_io_error_limit: G", str(cm.exception))
        self.assertEqual(self.counter.value, 5)
        self.assertEqual(self.sleeps, [])

    def test_node_limit_reached_refuses_to_read(self):
        self.counter.value = 5
        fp = self.plan([])
        with self.assertRaises(self.m.GenomeReadFail) as cm:
            self.read()
        self.assertIn("node_io_error_limit_reached", str(cm.exception))
        self.assertEqual(fp.attempts(self.p), 0)


class _NullLock:
    def __enter__(self):
        return self

    def __exit__(self, *a):
        pass


class TIOBuild(unittest.TestCase):
    """build() in-process (Pool workers are forked, so the patched scan_fasta and
    back-off are inherited) with faults injected per genome."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="t248iob.", dir=os.environ.get("TMPDIR"))
        self.f = Fixture(self.tmp)
        self.m = load_248()
        self.m.IO_BACKOFF_SECONDS = (0.0, 0.0, 0.0)
        self.old_env = os.environ.get("GIT_COMMIT")
        os.environ["GIT_COMMIT"] = self.f.commit

    def tearDown(self):
        if self.old_env is None:
            os.environ.pop("GIT_COMMIT", None)
        else:
            os.environ["GIT_COMMIT"] = self.old_env
        shutil.rmtree(self.tmp)

    def build(self, plan, out="run"):
        if plan is not None:
            fp = FaultPlan(os.path.join(self.tmp, "plan_" + out), plan)
            orig = self.m.scan_fasta
            self.m.scan_fasta = lambda path, label: orig(path, label, opener=fp.opener_for(path))
        r = self.f.root
        a = argparse.Namespace(
            queries_fna=os.path.join(r, "q.fna"), queries_tsv=os.path.join(r, "q.tsv"),
            corpus=os.path.join(r, "corpus.tsv"), expect=os.path.join(r, "expect.json"),
            blast_bin=self.f.blast, repo=self.f.repo, code_files=self.f.code_files,
            tooltests=self.f.tt, threads=2, out=os.path.join(r, out))
        with open(os.devnull, "w") as dn:
            old = sys.stdout, sys.stderr
            sys.stdout = sys.stderr = dn
            try:
                return self.m.build(a)
            finally:
                sys.stdout, sys.stderr = old

    def rd(self, out, name):
        with open(os.path.join(self.tmp, out, name), "rb") as f:
            return f.read()

    def test_transient_error_recovers_and_is_recorded(self):
        self.assertEqual(self.build(None, "clean"), 0)
        self.m = load_248()
        self.m.IO_BACKOFF_SECONDS = (0.0, 0.0, 0.0)
        self.assertEqual(self.build({"GCA_3.1.fna": [errno.EIO]}, "run"), 0)
        reads = self.rd("run", "genome_reads.tsv").decode().splitlines()
        self.assertEqual(reads[0], "shard\tassembly\trereads\tio_errors")
        self.assertIn("0\tGCA_3.1\t2\tEIO", reads)
        self.assertIn("1\tGCA_2.1\t0\t-", reads)
        res = json.loads(self.rd("run", "v0_result.json"))
        io = res["io_rereads"]
        self.assertEqual(io["node_io_errors"], 1)
        self.assertEqual(io["node"], res["host"])
        self.assertEqual(io["genomes_reread"], [{"assembly": "GCA_3.1", "rereads": 2, "io_errors": "EIO"}])
        self.assertEqual([(d["assembly"], d["errno"], d["attempt"], d["host"]) for d in io["events"]],
                         [("GCA_3.1", "EIO", 1, res["host"])])
        self.assertIn("genome_reads.tsv", res["outputs"])
        # the trusted manifest does not depend on I/O history
        for name in ("genomes.tsv", "contigs.tsv", "genome_manifest.index"):
            self.assertEqual(self.rd("run", name), self.rd("clean", name), name)
        self.assertEqual(self.f.check("run").returncode, 0)

    def test_persistent_error_fails_build(self):
        self.assertEqual(self.build({"GCA_3.1.fna": [errno.EIO] * 10}), 1)
        failed = self.rd("run", "FAILED").decode()
        self.assertIn("genome_io_retries_exhausted: GCA_3.1", failed)
        self.assertIn("host=", failed)
        self.assertFalse(os.path.exists(os.path.join(self.tmp, "run", "v0_result.json")))

    def test_node_limit_fails_build_naming_node(self):
        import socket
        plan = {"GCA_1.1.fna": [errno.EIO, errno.EIO], "GCA_2.1.fna": [errno.EIO, errno.EIO],
                "GCA_4.1.fna": [errno.ENOTCONN, errno.ENOTCONN]}
        self.assertEqual(self.build(plan), 1)
        failed = self.rd("run", "FAILED").decode()
        self.assertIn(f"node_io_error_limit: node {socket.gethostname()}", failed)
        self.assertIn(f"host={socket.gethostname()}", failed)
        # the record names the genomes and errnos even though the pool was torn down
        self.assertNotIn("logged errors (assembly:errno@time#attempt): []", failed)
        self.assertRegex(failed, r"GCA_[124]\.1:(EIO|ENOTCONN)@")
        self.assertIn("errno counts {", failed)
        self.assertFalse(os.path.exists(os.path.join(self.tmp, "run", "v0_result.json")))

    def test_reads_record_tamper_caught_by_check(self):
        self.assertEqual(self.build(None), 0)
        p = os.path.join(self.tmp, "run", "genome_reads.tsv")
        data = open(p).read().replace("0\tGCA_3.1\t0\t-", "0\tGCA_3.1\t1\t-")
        open(p, "w").write(data)
        r = self.f.check("run")
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("genome_reads_count: GCA_3.1", r.stderr)

    def test_consistent_reads_tamper_caught_by_result(self):
        # internally consistent edit (0/- -> 2/EIO): only v0_result.json can catch it
        self.assertEqual(self.build(None), 0)
        p = os.path.join(self.tmp, "run", "genome_reads.tsv")
        data = open(p).read().replace("0\tGCA_3.1\t0\t-", "0\tGCA_3.1\t2\tEIO")
        open(p, "w").write(data)
        r = self.f.check("run")
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("result_output_digest: genome_reads.tsv", r.stderr)

    def test_result_io_summary_tamper_caught(self):
        self.assertEqual(self.build(None), 0)
        p = os.path.join(self.tmp, "run", "v0_result.json")
        res = json.loads(open(p).read())
        res["io_rereads"]["genomes_reread"] = [{"assembly": "GCA_3.1", "rereads": 2, "io_errors": "EIO"}]
        open(p, "w").write(json.dumps(res))
        r = self.f.check("run")
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("result_io_rereads", r.stderr)

    def test_check_requires_result(self):
        self.assertEqual(self.build(None), 0)
        os.remove(os.path.join(self.tmp, "run", "v0_result.json"))
        r = self.f.check("run")
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("result_missing", r.stderr)


if __name__ == "__main__":
    unittest.main(verbosity=2)
