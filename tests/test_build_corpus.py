#!/usr/bin/env python3
"""Tests for scripts/250_build_corpus.py and tools/download_staging.sh (SPEC M2-C).

Hermetic: synthetic source directories, FASTA files and GenBank summary
in a temp dir; the download tool runs against a stub curl that
serves local files. Runs on the login node in seconds.
"""
import csv
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SCRIPT = os.path.join(ROOT, "scripts", "250_build_corpus.py")
DOWNLOAD = os.path.join(ROOT, "tools", "download_staging.sh")
COMMIT = "0" * 40

SEQ1 = "ACGTTGCAAGGCTTAACCGGTTAACGTACGATCGATCGTAGCTAGCTAGGATCC" * 3
SEQ1B = "TTGACCAGTAGGCATCGAT" * 2
SEQ2 = "GGGCCCATATATCGCGTTAAGGCCTTAAGG" * 4
SEQ3 = "CATCATGATGATCCGGAATTCCAAGGTTCC" * 4


def fasta(records, width=60, eol="\n", lower=False):
    out = []
    for name, seq in records:
        out.append(">" + name + eol)
        s = seq.lower() if lower else seq
        out.extend(s[i:i + width] + eol for i in range(0, len(s), width))
    return "".join(out)


def sha(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


SUMMARY_HDR = ("#assembly_accession\tbioproject\tspecies_taxid\tassembly_level\t"
               "seq_rel_date\tftp_path\tgbrs_paired_asm\tpaired_asm_comp")


def summary_line(acc, taxid, level="Complete Genome", date="2026-01-01",
                 paired="na", comp="na"):
    return "\t".join([acc, "PRJNA1", str(taxid), level, date,
                      "https://stub/%s" % acc, paired, comp])


class Tree:
    """Two source dirs (alpha: a GCA/GCF pair, a plain GCA, a plasmid-only
    sub-directory and a .fai; beta: a symlinked genome and another)."""

    def __init__(self, root):
        self.root = root
        self.a = os.path.join(root, "src_alpha")
        self.b = os.path.join(root, "src_beta")
        elsewhere = os.path.join(root, "elsewhere")
        for d in (self.a, self.b, elsewhere, os.path.join(self.a, "_plasmid_only")):
            os.makedirs(d)
        self.put(self.a, "GCA_000000001.1", fasta([("chr", SEQ1), ("p1", SEQ1B)]))
        # the RefSeq copy: other headers, record order, width, case, CRLF
        self.put(self.a, "GCF_000000001.1",
                 fasta([("NZ_p1 plasmid", SEQ1B), ("NZ_chr x", SEQ1)], width=17,
                       eol="\r\n", lower=True))
        self.put(self.a, "GCA_000000002.1", fasta([("c", SEQ2)]))
        self.put(os.path.join(self.a, "_plasmid_only"), "GCA_000000009.1", fasta([("p", SEQ1B)]))
        with open(os.path.join(self.a, "GCA_000000001.1.fna.fai"), "w") as fh:
            fh.write("chr\t162\t5\t60\t61\n")
        self.put(elsewhere, "GCA_000000003.1", fasta([("c", SEQ3)]))
        os.symlink(os.path.join(elsewhere, "GCA_000000003.1.fna"),
                   os.path.join(self.b, "GCA_000000003.1.fna"))
        self.put(self.b, "GCA_000000004.2", fasta([("c", SEQ2)]))      # same seq as ...002
        self.counts = {"alpha": 3, "beta": 2}
        self.write_sources()

    def put(self, d, acc, text):
        with open(os.path.join(d, acc + ".fna"), "w", newline="") as fh:
            fh.write(text)

    def write_sources(self, rows=None):
        self.sources = os.path.join(self.root, "sources.tsv")
        rows = rows or [["alpha", "101", "Alpha one", self.a, str(self.counts["alpha"])],
                        ["beta", "202", "Beta two", self.b, str(self.counts["beta"])]]
        with open(self.sources, "w") as fh:
            fh.write("species\tspecies_taxid\torganism\tsource_dir\texpected_n\n")
            fh.writelines("\t".join(r) + "\n" for r in rows)


def run(*args, commit=COMMIT, env_extra=None):
    env = dict(os.environ)
    env.pop("GIT_COMMIT", None)
    if commit:
        env["GIT_COMMIT"] = commit
    env.update(env_extra or {})
    return subprocess.run([sys.executable, SCRIPT] + list(args), capture_output=True,
                          text=True, env=env)


class Base(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="t250.", dir=os.environ.get("TMPDIR"))
        self.t = Tree(self.tmp)
        self.n = 0

    def tearDown(self):
        for dp, dn, fn in os.walk(self.tmp):
            for x in dn + fn:
                try:
                    os.chmod(os.path.join(dp, x), 0o755, follow_symlinks=False)
                except (OSError, NotImplementedError):
                    pass
        shutil.rmtree(self.tmp)

    def out(self):
        self.n += 1
        return os.path.join(self.tmp, "out%d" % self.n)

    def ok(self, r):
        self.assertEqual(r.returncode, 0, r.stderr)

    def fails(self, r, frag, out=None, outputs=()):
        self.assertEqual(r.returncode, 2, r.stderr)
        self.assertIn(frag, r.stderr)
        if out:
            for name in outputs:
                self.assertFalse(os.path.exists(os.path.join(out, name)), name)
            self.assertFalse(any(f.endswith(".json") and f != "FAILED.json"
                                 for f in os.listdir(out)))
            self.assertFalse(os.path.exists(os.path.join(out, ".tmp")))

    def build(self):
        o = self.out()
        self.ok(run("build", "--sources", self.t.sources, "--out", o))
        return os.path.join(o, "corpus.tsv")

    def rows(self, path):
        with open(path, newline="") as fh:
            return list(csv.DictReader(fh, delimiter="\t"))


class TBuild(Base):
    def test_rows_sorted_labelled_paths_unresolved(self):
        o = self.out()
        r = run("build", "--sources", self.t.sources, "--out", o)
        self.ok(r)
        rows = self.rows(os.path.join(o, "corpus.tsv"))
        self.assertEqual([x["assembly"] for x in rows],
                         ["GCA_000000001.1", "GCA_000000002.1", "GCA_000000003.1",
                          "GCA_000000004.2", "GCF_000000001.1"])
        self.assertEqual({x["assembly"]: x["species"] for x in rows}["GCA_000000003.1"], "beta")
        # the symlink path is kept as listed, not resolved; sub-dir and .fai ignored
        self.assertEqual(rows[2]["fna_path"], self.t.b + "/GCA_000000003.1.fna")
        self.assertNotIn("GCA_000000009.1", {x["assembly"] for x in rows})
        rec = json.load(open(os.path.join(o, "build.json")))
        self.assertEqual(rec["status"], "complete")
        self.assertEqual(rec["counts"]["assemblies"], 5)
        self.assertEqual(rec["counts"]["symlinks"], 1)
        self.assertEqual(rec["git_commit_submit"], COMMIT)
        self.assertEqual(rec["outputs"]["corpus.tsv"]["sha256"], sha(os.path.join(o, "corpus.tsv")))
        self.assertEqual(rec["inputs"]["sources"]["sha256"], sha(self.t.sources))
        raw = open(os.path.join(o, "corpus.tsv"), "rb").read()
        self.assertTrue(raw.startswith(b"assembly\tspecies\tfna_path\n") and raw.endswith(b"\n"))
        self.assertNotIn(b"\r", raw)

    def test_expect_sha256_pins_the_result(self):
        o = self.out()
        self.ok(run("build", "--sources", self.t.sources, "--out", o))
        good = sha(os.path.join(o, "corpus.tsv"))
        self.ok(run("build", "--sources", self.t.sources, "--expect-sha256", good,
                    "--out", self.out()))
        o3 = self.out()
        self.fails(run("build", "--sources", self.t.sources, "--expect-sha256", "f" * 64,
                       "--out", o3), "expected", o3, ["corpus.tsv"])
        self.assertTrue(os.path.exists(os.path.join(o3, "FAILED.json")))

    def test_count_must_equal_declared(self):
        self.t.counts["alpha"] = 4
        self.t.write_sources()
        o = self.out()
        self.fails(run("build", "--sources", self.t.sources, "--out", o),
                   "config declares 4", o, ["corpus.tsv"])

    def test_bad_entries_are_fatal(self):
        cases = {
            "broken symlink": lambda: os.symlink("/nonexistent/x.fna",
                                                 os.path.join(self.t.b, "GCA_000000005.1.fna")),
            "is empty": lambda: open(os.path.join(self.t.b, "GCA_000000005.1.fna"), "w").close(),
            "not an NCBI assembly accession": lambda: self.t.put(self.t.b, "genome5", ">c\nACGT\n"),
            "in two sources": lambda: self.t.put(self.t.b, "GCA_000000002.1", ">c\nACGT\n"),
            "is a directory": lambda: os.makedirs(os.path.join(self.t.b, "GCA_000000005.1.fna")),
        }
        for frag, make in cases.items():
            with self.subTest(frag):
                for x in os.listdir(self.t.b):
                    if x not in ("GCA_000000003.1.fna", "GCA_000000004.2.fna"):
                        p = os.path.join(self.t.b, x)
                        if os.path.isdir(p) and not os.path.islink(p):
                            os.rmdir(p)
                        else:
                            os.remove(p)
                make()
                o = self.out()
                self.fails(run("build", "--sources", self.t.sources, "--out", o), frag,
                           o, ["corpus.tsv"])

    def test_missing_source_dir_is_fatal(self):
        shutil.rmtree(self.t.b)
        o = self.out()
        self.fails(run("build", "--sources", self.t.sources, "--out", o), "does not exist",
                   o, ["corpus.tsv"])

    def test_bad_config_is_fatal(self):
        for rows, frag in (
                ([["alpha", "101", "A", self.t.a, "3"], ["alpha", "202", "B", self.t.b, "2"]],
                 "duplicate species"),
                ([["alpha", "101", "A", "relative/dir", "3"]], "invalid row"),
                ([["alpha", "x", "A", self.t.a, "3"]], "invalid row")):
            with self.subTest(frag):
                self.t.write_sources(rows)
                self.fails(run("build", "--sources", self.t.sources, "--out", self.out()), frag)

    def test_existing_outdir_refused_and_commit_required(self):
        o = self.out()
        os.makedirs(o)
        self.fails(run("build", "--sources", self.t.sources, "--out", o), "already exists")
        self.assertEqual(os.listdir(o), [])
        self.fails(run("build", "--sources", self.t.sources, "--out", self.out(), commit=None),
                   "GIT_COMMIT")


class TPairs(Base):
    def setUp(self):
        super().setUp()
        self.corpus = self.build()
        self.summary = os.path.join(self.tmp, "summary.txt")
        self.write_summary([summary_line("GCA_000000001.1", 101, paired="GCF_000000001.1",
                                         comp="identical")])

    def write_summary(self, lines, header=SUMMARY_HDR):
        with open(self.summary, "w") as fh:
            fh.write("##  comment line\n" + header + "\n")
            fh.writelines(l + "\n" for l in lines)

    def pairs(self, *extra, out=None):
        o = out or self.out()
        return o, run("pairs", "--corpus", self.corpus, "--corpus-sha256", sha(self.corpus),
                      "--summary", self.summary, "--summary-sha256", sha(self.summary),
                      "--workers", "2", "--out", o, *extra)

    def test_identical_pair_collapses_into_the_gca(self):
        o, r = self.pairs()
        self.ok(r)
        d = self.rows(os.path.join(o, "corpus_duplicates.tsv"))
        self.assertEqual(len(d), 1)
        self.assertEqual((d[0]["duplicate"], d[0]["canonical"], d[0]["relation"],
                          d[0]["collapse"]),
                         ("GCF_000000001.1", "GCA_000000001.1", "identical_sequence", "yes"))
        self.assertEqual((d[0]["ncbi_paired_asm"], d[0]["ncbi_paired_comp"]),
                         ("GCF_000000001.1", "identical"))
        self.assertEqual(d[0]["dup_bases"], str(len(SEQ1) + len(SEQ1B)))
        rec = json.load(open(os.path.join(o, "pairs.json")))
        self.assertEqual(rec["counts"]["assemblies_after_collapse"], 4)
        self.assertEqual(rec["counts"]["genomes_read"], 5)
        self.assertEqual(rec["corpus_sha256"], sha(self.corpus))
        for k in ("corpus", "summary"):
            self.assertIn("sha256", rec["inputs"][k])

    def test_every_genome_read_is_recorded(self):
        o, r = self.pairs()
        self.ok(r)
        content = self.rows(os.path.join(o, "corpus_content.tsv"))
        corpus = self.rows(self.corpus)
        self.assertEqual([c["assembly"] for c in content], [c["assembly"] for c in corpus])
        for c in content:
            self.assertEqual(c["file_sha256"], sha(c["fna_path"]))
            self.assertEqual(int(c["bytes"]), os.path.getsize(c["fna_path"]))
            self.assertEqual(c["realpath"], os.path.realpath(c["fna_path"]))
        sym = {c["assembly"]: c for c in content}["GCA_000000003.1"]
        self.assertNotEqual(sym["realpath"], sym["fna_path"])

    def test_every_ascii_whitespace_is_ignored(self):
        p = os.path.join(self.t.a, "GCF_000000001.1.fna")
        text = open(p, newline="").read().replace("\r\n", "\x0b\r\n", 3).replace("a", "a\x0c", 1)
        open(p, "w", newline="").write(text)
        self.corpus = self.build()
        o, r = self.pairs()
        self.ok(r)
        d = self.rows(os.path.join(o, "corpus_duplicates.tsv"))[0]
        self.assertEqual(d["relation"], "identical_sequence")

    def test_one_base_difference_is_not_collapsed(self):
        p = os.path.join(self.t.a, "GCF_000000001.1.fna")
        text = open(p, newline="").read()
        i = text.index("\r\n") + 2                 # first base of the first record
        text = text[:i] + ("c" if text[i] != "c" else "g") + text[i + 1:]
        open(p, "w", newline="").write(text)
        self.corpus = self.build()
        o, r = self.pairs()
        self.ok(r)
        d = self.rows(os.path.join(o, "corpus_duplicates.tsv"))[0]
        self.assertEqual((d["relation"], d["collapse"]), ("differs", "no"))

    def test_other_version_collapses_as_different_version(self):
        os.rename(os.path.join(self.t.a, "GCA_000000001.1.fna"),
                  os.path.join(self.t.a, "GCA_000000001.2.fna"))
        os.remove(os.path.join(self.t.a, "GCA_000000001.1.fna.fai"))
        self.t.put(self.t.a, "GCA_000000001.2", fasta([("chr", SEQ1 + "A")]))
        self.corpus = self.build()
        o, r = self.pairs()
        self.ok(r)
        d = self.rows(os.path.join(o, "corpus_duplicates.tsv"))[0]
        self.assertEqual((d["canonical"], d["relation"], d["collapse"]),
                         ("GCA_000000001.2", "different_version", "yes"))
        self.assertEqual(d["ncbi_paired_asm"], "not_in_summary")

    def test_identical_sequence_under_other_ids_is_listed_not_collapsed(self):
        o, r = self.pairs()
        self.ok(r)
        ident = self.rows(os.path.join(o, "corpus_identical_other_id.tsv"))
        self.assertEqual({x["assembly"] for x in ident}, {"GCA_000000002.1", "GCA_000000004.2"})
        self.assertEqual({x["group_id"] for x in ident}, {"I00001"})
        dups = {x["duplicate"] for x in self.rows(os.path.join(o, "corpus_duplicates.tsv"))}
        self.assertFalse(dups & {"GCA_000000002.1", "GCA_000000004.2"})
        # same contig lengths but one base different: not a group
        self.t.put(self.t.b, "GCA_000000004.2", fasta([("c", "T" + SEQ2[1:])]))
        self.corpus = self.build()
        o, r = self.pairs()
        self.ok(r)
        self.assertEqual(self.rows(os.path.join(o, "corpus_identical_other_id.tsv")), [])

    def test_two_genbank_versions_keep_the_newest(self):
        self.t.put(self.t.b, "GCA_000000006.1", fasta([("c", SEQ3 + "AC")]))
        self.t.put(self.t.b, "GCA_000000006.2", fasta([("c", SEQ3 + "GT")]))
        self.t.counts["beta"] = 4
        self.t.write_sources()
        self.corpus = self.build()
        o, r = self.pairs()
        self.ok(r)
        d = {x["duplicate"]: x for x in self.rows(os.path.join(o, "corpus_duplicates.tsv"))}
        self.assertEqual((d["GCA_000000006.1"]["canonical"], d["GCA_000000006.1"]["relation"]),
                         ("GCA_000000006.2", "different_version"))

    def test_same_id_pair_is_not_an_other_id_group(self):
        # GCA_1 and GCF_1 are identical but ONE id: never an other-id group;
        # a third id with the same sequence joins them into one
        o, r = self.pairs()
        self.ok(r)
        ident = self.rows(os.path.join(o, "corpus_identical_other_id.tsv"))
        self.assertEqual({x["assembly"] for x in ident}, {"GCA_000000002.1", "GCA_000000004.2"})
        self.t.put(self.t.b, "GCA_000000005.1", fasta([("chr", SEQ1), ("p", SEQ1B)]))
        self.t.counts["beta"] = 3
        self.t.write_sources()
        self.corpus = self.build()
        o, r = self.pairs()
        self.ok(r)
        ident = self.rows(os.path.join(o, "corpus_identical_other_id.tsv"))
        g = {x["assembly"]: x for x in ident if x["assembly"].endswith("01.1")
             or x["assembly"] == "GCA_000000005.1"}
        self.assertEqual(set(g), {"GCA_000000001.1", "GCF_000000001.1", "GCA_000000005.1"})
        self.assertEqual({x["distinct_ids"] for x in g.values()}, {"2"})

    def test_pinned_inputs(self):
        o = self.out()
        r = run("pairs", "--corpus", self.corpus, "--corpus-sha256", "a" * 64,
                "--summary", self.summary, "--summary-sha256", sha(self.summary), "--out", o)
        self.fails(r, "expected", o, ["corpus_duplicates.tsv"])

    def test_unreadable_genome_fails_without_outputs(self):
        os.chmod(os.path.join(self.t.a, "GCF_000000001.1.fna"), 0)
        if os.access(os.path.join(self.t.a, "GCF_000000001.1.fna"), os.R_OK):
            self.skipTest("running as a user who can read mode-000 files")
        o, r = self.pairs()
        self.fails(r, "cannot read GCF_000000001.1", o,
                   ["corpus_duplicates.tsv", "corpus_identical_other_id.tsv"])

    def test_git_state_must_be_measurable(self):
        o = self.out()
        r = run("pairs", "--corpus", self.corpus, "--corpus-sha256", sha(self.corpus),
                "--summary", self.summary, "--summary-sha256", sha(self.summary), "--out", o,
                env_extra={"PATH": "/nonexistent"})       # no git on PATH
        self.fails(r, "git_commit_runtime")
        self.assertFalse(os.path.exists(o))

    def test_malformed_summary_line_is_fatal(self):
        self.write_summary([summary_line("GCA_000000001.1", 101) + "\textra"])
        o, r = self.pairs()
        self.fails(r, "fields, header has", o, ["corpus_duplicates.tsv"])


class TInProcess(Base):
    """Failures injected inside 250 (short read, failing publish)."""

    def setUp(self):
        super().setUp()
        import importlib.util
        spec = importlib.util.spec_from_file_location("m250", SCRIPT)
        self.m = importlib.util.module_from_spec(spec)
        sys.modules["m250"] = self.m            # the worker pool pickles by module name
        spec.loader.exec_module(self.m)
        self.corpus = self.build()
        self.summary = os.path.join(self.tmp, "summary.txt")
        with open(self.summary, "w") as fh:
            fh.write(SUMMARY_HDR + "\n" + summary_line("GCA_000000001.1", 101) + "\n")
        os.environ["GIT_COMMIT"] = COMMIT

    def tearDown(self):
        os.environ.pop("GIT_COMMIT", None)
        super().tearDown()

    def argv(self, o):
        return ["pairs", "--corpus", self.corpus, "--corpus-sha256", sha(self.corpus),
                "--summary", self.summary, "--summary-sha256", sha(self.summary),
                "--workers", "1", "--out", o]

    def test_short_read_is_fatal(self):
        p = self.rows(self.corpus)[0]["fna_path"]
        real = self.m.os.fstat

        class Big:
            def __init__(self, st):
                self.st_size = st.st_size + 1        # the file claims one byte more
        self.m.os.fstat = lambda fd: Big(real(fd))
        try:
            with self.assertRaises(self.m.Fail) as cm:
                self.m.read_genome(p)
        finally:
            self.m.os.fstat = real
        self.assertIn("short read", str(cm.exception))

    def test_failure_while_publishing_leaves_nothing(self):
        o = self.out()
        real, calls = os.replace, []

        def flaky(src, dst):
            calls.append(dst)
            if len(calls) == 2:                     # the second output's move
                raise OSError(5, "Input/output error")
            return real(src, dst)
        self.m.os.replace = flaky
        try:
            rc = self.m.main(self.argv(o))
        finally:
            self.m.os.replace = real
        self.assertEqual(rc, 2)
        self.assertGreaterEqual(len(calls), 2)
        self.assertEqual(sorted(os.listdir(o)), ["FAILED.json"])

    def test_input_short_read_is_fatal(self):
        # [M2R4-1] the input reader checks bytes read against the file size
        real = self.m.os.fstat

        class Big:
            def __init__(self, st):
                self.st_size = st.st_size + 1
        self.m.os.fstat = lambda fd: Big(real(fd))
        try:
            with self.assertRaises(self.m.Fail) as cm:
                list(self.m.Input(self.corpus, "corpus").lines())
        finally:
            self.m.os.fstat = real
        self.assertIn("short read", str(cm.exception))

    def opened_once(self, argv, tables):
        """Run main(argv) with every second (or later) open of an input table
        served a copy cut at a line boundary -- the parse must never see it."""
        counts = {p: 0 for p in tables}
        cut = {}
        for p in tables:
            data = open(p, "rb").read()
            cut[p] = os.path.join(self.tmp, "cut_" + str(len(cut)))
            with open(cut[p], "wb") as fh:
                fh.write(data[:data.rindex(b"\n", 0, len(data) - 1) + 1])
        real_open = open

        def spy(path, *a, **k):
            if path in counts:
                counts[path] += 1
                if counts[path] > 1:
                    path = cut[path]
            return real_open(path, *a, **k)
        self.m.open = spy
        try:
            rc = self.m.main(argv)
        finally:
            del self.m.open
        return rc, counts

    def test_each_input_table_is_read_once(self):
        o = self.out()
        rc, counts = self.opened_once(self.argv(o), [self.corpus, self.summary])
        self.assertEqual(rc, 0)
        self.assertEqual(counts, {self.corpus: 1, self.summary: 1})
        rec = json.load(open(os.path.join(o, "pairs.json")))
        self.assertEqual(rec["inputs"]["corpus"]["sha256"], sha(self.corpus))
        self.assertEqual(rec["inputs"]["summary"]["sha256"], sha(self.summary))
        self.assertEqual(len(self.rows(os.path.join(o, "corpus_content.tsv"))),
                         len(self.rows(self.corpus)))
        o = self.out()
        rc, counts = self.opened_once(
            ["gaps", "--sources", self.t.sources, "--corpus", self.corpus,
             "--corpus-sha256", sha(self.corpus), "--summary", self.summary,
             "--summary-sha256", sha(self.summary), "--cutoff", "2026-09-01", "--out", o],
            [self.t.sources, self.corpus, self.summary])
        self.assertEqual(rc, 0)
        self.assertEqual(set(counts.values()), {1})
        rec = json.load(open(os.path.join(o, "gaps.json")))
        self.assertEqual(rec["inputs"]["sources"]["sha256"], sha(self.t.sources))
        o = self.out()
        rc, counts = self.opened_once(["build", "--sources", self.t.sources, "--out", o],
                                      [self.t.sources])
        self.assertEqual((rc, counts), (0, {self.t.sources: 1}))

    def test_failure_creating_tmp_leaves_failed_json(self):
        # [M2R3-5] the run directory exists, .tmp cannot be created
        o = self.out()
        real = os.makedirs

        def no_tmp(path, *a, **k):
            if os.path.basename(path) == ".tmp":
                raise OSError(28, "No space left on device")
            return real(path, *a, **k)
        self.m.os.makedirs = no_tmp
        try:
            rc = self.m.main(self.argv(o))
        finally:
            self.m.os.makedirs = real
        self.assertEqual(rc, 2)
        self.assertEqual(sorted(os.listdir(o)), ["FAILED.json"])
        rec = json.load(open(os.path.join(o, "FAILED.json")))
        self.assertEqual(rec["status"], "failed")
        self.assertIn("No space left", rec["reason"])

    def test_refused_directory_is_never_touched(self):
        # a directory this run did not create gets no FAILED.json either,
        # including one created between the check and makedirs
        o = self.out()
        real = os.makedirs

        def racer(path, *a, **k):
            if path == o:
                real(path)                          # someone else got there first
            return real(path, *a, **k)
        self.m.os.makedirs = racer
        try:
            rc = self.m.main(self.argv(o))
        finally:
            self.m.os.makedirs = real
        self.assertEqual(rc, 2)
        self.assertEqual(os.listdir(o), [])


class TGaps(Base):
    def setUp(self):
        super().setUp()
        self.corpus = self.build()
        self.summary = os.path.join(self.tmp, "summary.txt")

    def gaps(self, lines, cutoff="2026-09-01", header=SUMMARY_HDR, pre=""):
        with open(self.summary, "w") as fh:
            fh.write(pre + "##  comment\n" + header + "\n")
            fh.writelines(l + "\n" for l in lines)
        o = self.out()
        r = run("gaps", "--sources", self.t.sources, "--corpus", self.corpus,
                "--corpus-sha256", sha(self.corpus), "--summary", self.summary,
                "--summary-sha256", sha(self.summary), "--cutoff", cutoff, "--out", o)
        return o, r

    def test_statuses(self):
        lines = [
            summary_line("GCA_000000001.1", 101),                         # in corpus
            summary_line("GCA_000000002.3", 101),                         # newer version
            summary_line("GCA_000000009.1", 101),                         # plasmid-only aside
            summary_line("GCA_000000010.1", 101, date="2026-09-01"),      # missing (= cutoff)
            summary_line("GCA_000000011.1", 101, date="2026-09-02"),      # released after
            summary_line("GCA_000000012.1", 101, level="Contig"),         # not usable
            summary_line("GCA_000000013.1", 999),                         # other species
            summary_line("GCA_000000003.1", 202),
            # beta's GCA_000000004.2 is absent: not in the usable set any more
        ]
        o, r = self.gaps(lines)
        self.ok(r)
        per = {x["species"]: x for x in self.rows(os.path.join(o, "corpus_gaps.tsv"))}
        a, b = per["alpha"], per["beta"]
        self.assertEqual((a["corpus_assemblies"], a["corpus_gcf"], a["plasmid_only_set_aside"],
                          a["usable_total"], a["usable_at_cutoff"], a["missing_at_cutoff"],
                          a["released_after_cutoff"], a["version_mismatch"],
                          a["not_in_usable_set"]),
                         ("3", "1", "1", "5", "4", "1", "1", "1", "0"))
        self.assertEqual((b["usable_total"], b["not_in_usable_set"], b["missing_at_cutoff"]),
                         ("1", "1", "0"))
        gl = {(x["accession"], x["status"]) for x in self.rows(os.path.join(o, "corpus_gap_list.tsv"))}
        self.assertEqual(gl, {("GCA_000000002.3", "version_mismatch"),
                              ("GCA_000000010.1", "missing_at_cutoff"),
                              ("GCA_000000011.1", "released_after_cutoff"),
                              ("GCA_000000004.2", "not_in_usable_set")})
        rec = json.load(open(os.path.join(o, "gaps.json")))
        self.assertEqual(rec["cutoff"], "2026-09-01")
        self.assertEqual(rec["counts"]["summary_data_lines"], 8)

    def test_queue_holds_only_downloadable_missing_rows(self):
        lines = [summary_line("GCA_000000010.1", 101),
                 summary_line("GCA_000000011.1", 101, date="2026-10-01"),
                 summary_line("GCA_000000002.3", 101),
                 "\t".join(["GCA_000000012.1", "PRJNA1", "101", "Chromosome", "2026-01-01",
                            "na", "na", "na"]),
                 # [M2R3-2] no host
                 "\t".join(["GCA_000000014.1", "PRJNA1", "101", "Chromosome", "2026-01-01",
                            "https:///GCA_000000014.1", "na", "na"]),
                 summary_line("GCA_000000003.1", 202)]
        o, r = self.gaps(lines)
        self.ok(r)
        q = self.rows(os.path.join(o, "corpus_download_queue.tsv"))
        self.assertEqual([x["accession"] for x in q], ["GCA_000000010.1"])
        rec = json.load(open(os.path.join(o, "gaps.json")))
        self.assertEqual((rec["counts"]["queued"], rec["counts"]["missing_without_ftp_path"]), (1, 2))
        gl = {x["accession"] for x in self.rows(os.path.join(o, "corpus_gap_list.tsv"))}
        self.assertLessEqual({"GCA_000000012.1", "GCA_000000014.1"}, gl)   # reported, not queued

    def test_plasmid_only_listing_is_published(self):
        o, r = self.gaps([summary_line("GCA_000000009.1", 101)])
        self.ok(r)
        pl = self.rows(os.path.join(o, "corpus_plasmid_only.tsv"))
        self.assertEqual([(x["species"], x["accession"]) for x in pl], [("alpha", "GCA_000000009.1")])
        rec = json.load(open(os.path.join(o, "gaps.json")))
        self.assertEqual(rec["outputs"]["corpus_plasmid_only.tsv"]["sha256"],
                         sha(os.path.join(o, "corpus_plasmid_only.tsv")))

    def test_set_aside_versions_are_all_counted(self):
        # [M2R3-4] two versions of one id in _plasmid_only: both listed, both
        # counted; the summary's version among them is no mismatch
        self.t.put(os.path.join(self.t.a, "_plasmid_only"), "GCA_000000009.2", fasta([("p", SEQ1B)]))
        for summ_acc, mismatch in (("GCA_000000009.2", None),
                                   ("GCA_000000009.3", "GCA_000000009.1,GCA_000000009.2")):
            with self.subTest(summ_acc):
                o, r = self.gaps([summary_line(summ_acc, 101)])
                self.ok(r)
                pl = self.rows(os.path.join(o, "corpus_plasmid_only.tsv"))
                self.assertEqual([x["accession"] for x in pl], ["GCA_000000009.1", "GCA_000000009.2"])
                per = {x["species"]: x for x in self.rows(os.path.join(o, "corpus_gaps.tsv"))}
                self.assertEqual(per["alpha"]["plasmid_only_set_aside"], "2")
                rec = json.load(open(os.path.join(o, "gaps.json")))
                self.assertEqual(rec["counts"]["plasmid_only_listed"], 2)
                self.assertEqual(rec["counts"]["totals"]["plasmid_only_set_aside"], 2)
                gl = [x for x in self.rows(os.path.join(o, "corpus_gap_list.tsv"))
                      if x["status"] == "version_mismatch"]
                self.assertEqual([x["corpus_accession"] for x in gl], [mismatch] if mismatch else [])

    def test_uncheckable_plasmid_dir_is_fatal(self):
        pd = os.path.join(self.t.a, "_plasmid_only")
        os.chmod(pd, 0)
        try:
            if os.access(pd, os.R_OK):
                self.skipTest("running as a user who can read mode-000 directories")
            o, r = self.gaps([summary_line("GCA_000000010.1", 101)])
            self.fails(r, "_plasmid_only", o, ["corpus_gaps.tsv", "corpus_download_queue.tsv"])
        finally:
            os.chmod(pd, 0o755)
        shutil.rmtree(pd)
        open(pd, "w").close()                   # a file where the directory should be
        o, r = self.gaps([summary_line("GCA_000000010.1", 101)])
        self.fails(r, "not a directory", o, ["corpus_gaps.tsv"])
        os.remove(pd)
        os.symlink(os.path.join(self.tmp, "gone"), pd)      # [M2R3-3] a broken symlink
        o, r = self.gaps([summary_line("GCA_000000010.1", 101)])
        self.fails(r, "exists but cannot be checked", o,
                   ["corpus_gaps.tsv", "corpus_plasmid_only.tsv"])

    def test_bad_summary_or_cutoff_is_fatal(self):
        for kw, frag in (
                (dict(lines=[summary_line("GCA_000000010.1", 101)[:-3]]), "fields, header has"),
                (dict(lines=[summary_line("GCA_000000010.1", 101, date="2026/01/01")]),
                 "expected YYYY-MM-DD"),
                (dict(lines=[summary_line("GCF_000000010.1", 101)]), "not a GenBank"),
                (dict(lines=[summary_line("GCA_000000010.1", 101)],
                      pre=summary_line("GCA_000000020.1", 101) + "\n"), "data before"),
                (dict(lines=[summary_line("GCA_000000010.1", 101)], cutoff="2026/09/01"),
                 "--cutoff must be"),
                (dict(lines=[summary_line("GCA_000000010.1", 101)],
                      header=SUMMARY_HDR.replace("\tgbrs_paired_asm", "")), "header lacks")):
            with self.subTest(frag):
                o, r = self.gaps(**kw)
                self.fails(r, frag, o, ["corpus_gaps.tsv", "corpus_gap_list.tsv"])

    def test_corpus_species_must_be_declared(self):
        self.t.write_sources([["alpha", "101", "A", self.t.a, "3"]])
        o, r = self.gaps([summary_line("GCA_000000010.1", 101)])
        self.fails(r, "not declared", o, ["corpus_gaps.tsv"])


CURL_STUB = r'''#!/usr/bin/env python3
import os, sys
args, out, url = sys.argv[1:], None, None
for i, a in enumerate(args):
    if a == "-o":
        out = args[i + 1]
    elif a.startswith("http"):
        url = a
with open(os.environ["CURL_LOG"], "a") as f:
    f.write(url + "\n")
path = os.path.join(os.environ["CURL_ROOT"], url.split("://", 1)[1])
if not os.path.exists(path):
    sys.exit(22)
data = open(path, "rb").read()
if out:
    open(out, "wb").write(data)
else:
    sys.stdout.buffer.write(data)
if os.environ.get("CURL_BROKEN") and os.environ["CURL_BROKEN"] in url:
    sys.exit(18)                    # everything written, then a transfer error
'''


class TDownload(unittest.TestCase):
    """tools/download_staging.sh against a stub curl serving local files."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="tdl.", dir=os.environ.get("TMPDIR"))
        self.serve = os.path.join(self.tmp, "serve")
        self.bin = os.path.join(self.tmp, "bin")
        os.makedirs(self.bin)
        with open(os.path.join(self.bin, "curl"), "w") as fh:
            fh.write(CURL_STUB.replace("/usr/bin/env python3", sys.executable, 1))
        os.chmod(os.path.join(self.bin, "curl"), 0o755)
        self.log = os.path.join(self.tmp, "curl.log")
        self.queue = os.path.join(self.tmp, "queue.tsv")
        self.stage = os.path.join(self.tmp, "stage")

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def publish(self, acc, text, bad_md5=False):
        import gzip
        base = acc + "_ASM1v1"
        d = os.path.join(self.serve, "stub", "all", base)
        os.makedirs(d)
        gz = os.path.join(d, base + "_genomic.fna.gz")
        with gzip.open(gz, "wt") as fh:
            fh.write(text)
        md5 = hashlib.md5(open(gz, "rb").read()).hexdigest()
        if bad_md5:
            md5 = "0" * 32
        with open(os.path.join(d, "md5checksums.txt"), "w") as fh:
            fh.write("%s  ./%s_genomic.fna.gz\n" % (md5, base))
        return "https://stub/all/" + base + "/"

    def run_dl(self, rows, **env_extra):
        with open(self.queue, "w") as fh:
            fh.write("species\taccession\tstatus\tseq_rel_date\tassembly_level\tftp_path\n")
            fh.writelines("%s\t%s\tmissing_at_cutoff\t2026-01-01\tComplete Genome\t%s\n" % r
                          for r in rows)
        env = dict(os.environ, PATH=self.bin + os.pathsep + os.environ["PATH"],
                   CURL_LOG=self.log, CURL_ROOT=self.serve, TRIES="2", BACKOFF="0", PAR="2")
        env.update(env_extra)
        return subprocess.run(["bash", DOWNLOAD, self.queue, self.stage],
                              capture_output=True, text=True, env=env)

    def test_verified_download_then_resume_skips(self):
        url = self.publish("GCA_000000001.1", ">c\nACGT\n")
        r = self.run_dl([("alpha", "GCA_000000001.1", url)])
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        fna = os.path.join(self.stage, "alpha", "GCA_000000001.1.fna")
        self.assertEqual(open(fna).read(), ">c\nACGT\n")
        ok = dict(l.rstrip("\n").split("\t") for l in open(fna + ".ok"))
        self.assertEqual(ok["fna_sha256"], sha(fna))
        self.assertNotIn("//GCA", ok["url"])
        n = len(open(self.log).read().splitlines())
        r = self.run_dl([("alpha", "GCA_000000001.1", url)])
        self.assertEqual(r.returncode, 0)
        self.assertEqual(len(open(self.log).read().splitlines()), n)    # nothing re-fetched

    def test_incomplete_ok_record_is_not_accepted(self):
        # [M2R4-2] every field of the .ok record is checked, bytes against the .fna
        url = self.publish("GCA_000000001.1", ">c\nACGT\n")
        row = [("alpha", "GCA_000000001.1", url)]
        self.assertEqual(self.run_dl(row).returncode, 0)
        fna = os.path.join(self.stage, "alpha", "GCA_000000001.1.fna")
        good = open(fna + ".ok").read()
        lines = good.splitlines(True)
        tampered = {
            "no utc": "".join(l for l in lines if not l.startswith("utc")),
            "no fna_bytes": "".join(l for l in lines if not l.startswith("fna_bytes")),
            "bytes differ": good.replace("fna_bytes\t8\n", "fna_bytes\t9\n"),
            "other accession": good.replace("accession\tGCA_000000001.1", "accession\tGCA_000000002.1"),
            "bad utc": good.replace("Z\n", "\n"),
            "line twice": good + lines[2],
            "extra line": good + "note\tx\n"}
        self.assertNotEqual(tampered["bytes differ"], good)
        for name, text in tampered.items():
            with self.subTest(name):
                with open(fna + ".ok", "w") as fh:
                    fh.write(text)
                n = len(open(self.log).read().splitlines())
                r = self.run_dl(row)
                self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
                self.assertGreater(len(open(self.log).read().splitlines()), n)   # fetched again
                self.assertEqual(open(fna + ".ok").read().count("\n"), 6)
        # and the final check: an incomplete record with no server is MISSING
        with open(fna + ".ok", "w") as fh:
            fh.write(tampered["no utc"])
        shutil.rmtree(self.serve)
        r = self.run_dl(row)
        self.assertEqual(r.returncode, 1, r.stdout)
        self.assertIn("MISSING alpha GCA_000000001.1", r.stdout)

    def test_unmeasurable_counts_are_never_complete(self):
        # [M2R4-3] wc failing on download.fail (or on the checked queue)
        url = self.publish("GCA_000000001.1", ">c\nACGT\n")
        for target, rc, frag in (("download.fail", 1, "cannot read"),
                                 (".queue.checked", 2, "rejected")):
            with self.subTest(target):
                with open(os.path.join(self.bin, "wc"), "w") as fh:
                    fh.write('#!/bin/sh\ncase "$(readlink /proc/$$/fd/0)" in *%s) exit 1;; esac\n'
                             'exec /usr/bin/wc "$@"\n' % target)
                os.chmod(os.path.join(self.bin, "wc"), 0o755)
                r = self.run_dl([("alpha", "GCA_000000001.1", url)])
                self.assertEqual(r.returncode, rc, r.stdout + r.stderr)
                self.assertIn(frag, r.stdout)
                self.assertNotIn("COMPLETE\n", r.stdout.replace("INCOMPLETE", ""))

    def test_bad_settings_are_refused(self):
        url = self.publish("GCA_000000001.1", ">c\nACGT\n")
        for kw in ({"PAR": "0"}, {"TRIES": "x"}, {"BACKOFF": "-1"}):
            with self.subTest(kw):
                r = self.run_dl([("alpha", "GCA_000000001.1", url)], **kw)
                self.assertEqual(r.returncode, 2, r.stdout)
                self.assertFalse(os.path.exists(self.log))

    def test_md5_mismatch_is_rejected(self):
        good = self.publish("GCA_000000001.1", ">c\nACGT\n")
        bad = self.publish("GCA_000000002.1", ">c\nGGGG\n", bad_md5=True)
        r = self.run_dl([("alpha", "GCA_000000001.1", good), ("alpha", "GCA_000000002.1", bad)])
        self.assertEqual(r.returncode, 1, r.stdout)
        self.assertIn("INCOMPLETE", r.stdout)
        self.assertFalse(os.path.exists(os.path.join(self.stage, "alpha", "GCA_000000002.1.fna")))
        self.assertIn("FAIL alpha GCA_000000002.1", open(os.path.join(self.stage, "download.fail")).read())
        self.assertTrue(os.path.exists(os.path.join(self.stage, "alpha", "GCA_000000001.1.fna.ok")))


    def test_failed_md5_list_fetch_is_a_failure(self):
        # [M2R3-1] curl delivers the right md5 line, then exits non-zero
        url = self.publish("GCA_000000001.1", ">c\nACGT\n")
        r = self.run_dl([("alpha", "GCA_000000001.1", url)], CURL_BROKEN="md5checksums.txt")
        self.assertEqual(r.returncode, 1, r.stdout + r.stderr)
        self.assertIn("MISSING alpha GCA_000000001.1", r.stdout)
        self.assertNotIn("_genomic.fna.gz", open(self.log).read())     # never got that far
        self.assertEqual(os.listdir(os.path.join(self.stage, "alpha")), [])
        self.assertIn("FAIL alpha GCA_000000001.1", open(os.path.join(self.stage, "download.fail")).read())

    def test_two_md5_lines_for_the_file_are_a_failure(self):
        url = self.publish("GCA_000000001.1", ">c\nACGT\n")
        sums = os.path.join(self.serve, "stub", "all", "GCA_000000001.1_ASM1v1", "md5checksums.txt")
        line = open(sums).read()
        with open(sums, "a") as fh:
            fh.write(line)
        r = self.run_dl([("alpha", "GCA_000000001.1", url)])
        self.assertEqual(r.returncode, 1, r.stdout + r.stderr)
        self.assertNotIn("_genomic.fna.gz", open(self.log).read())

    def test_url_rule_matches_250(self):
        # [M2R3-2] one URL rule in both tools: a host name is required
        import importlib.util
        spec = importlib.util.spec_from_file_location("m250_url", SCRIPT)
        m = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(m)
        a = "GCA_000000001.1"
        good = ["https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/000/001/%s_ASM1v1" % a,
                "http://stub/%s" % a, "https://stub:8080/x/%s/" % a, "https://a.b-c9.d/%s" % a]
        bad = ["https:///%s" % a, "https://:80/%s" % a, "https://-stub/%s" % a,
               "https://stub-/%s" % a, "https://stub:/%s" % a, "https://user@stub/%s" % a,
               "https://%s" % a, "ftp://stub/%s" % a, "https://stub/%s?x=1" % a,
               "https://st_ub/%s" % a]
        for url in good + bad:
            with self.subTest(url):
                shutil.rmtree(self.stage, ignore_errors=True)
                r = self.run_dl([("alpha", a, url)], TRIES="1")
                accepted_250 = bool(m.FTP.fullmatch(url))
                self.assertEqual(accepted_250, url in good)
                self.assertEqual(r.returncode != 2, accepted_250, r.stdout + r.stderr)
                if not accepted_250:
                    self.assertIn("ftp_path", r.stderr)

    def test_stray_ok_outside_the_queue_does_not_mask_a_failure(self):
        stray = os.path.join(self.stage, "alpha")
        os.makedirs(stray)
        with open(os.path.join(stray, "GCA_000000099.1.fna"), "w") as fh:
            fh.write(">c\nA\n")
        with open(os.path.join(stray, "GCA_000000099.1.fna.ok"), "w") as fh:
            fh.write("gz_md5\t%s\nfna_sha256\t%s\n" % ("0" * 32, "0" * 64))
        good = self.publish("GCA_000000001.1", ">c\nACGT\n")
        bad = self.publish("GCA_000000002.1", ">c\nGGGG\n", bad_md5=True)
        r = self.run_dl([("alpha", "GCA_000000001.1", good), ("alpha", "GCA_000000002.1", bad)])
        self.assertEqual(r.returncode, 1, r.stdout)
        self.assertIn("MISSING alpha GCA_000000002.1", r.stdout)
        self.assertIn("verified 1", r.stdout)

    def test_failed_digest_writes_no_record(self):
        with open(os.path.join(self.bin, "sha256sum"), "w") as fh:
            fh.write("#!/bin/sh\nexit 1\n")
        os.chmod(os.path.join(self.bin, "sha256sum"), 0o755)
        url = self.publish("GCA_000000001.1", ">c\nACGT\n")
        r = self.run_dl([("alpha", "GCA_000000001.1", url)])
        self.assertEqual(r.returncode, 1, r.stdout)
        d = os.path.join(self.stage, "alpha")
        self.assertFalse(os.path.exists(os.path.join(d, "GCA_000000001.1.fna.ok")))
        self.assertFalse(os.path.exists(os.path.join(d, "GCA_000000001.1.fna")))

    def test_queue_header_is_checked(self):
        with open(self.queue, "w") as fh:
            fh.write("acc\turl\nGCA_000000001.1\thttps://stub/x\n")
        env = dict(os.environ, PATH=self.bin + os.pathsep + os.environ["PATH"],
                   CURL_LOG=self.log, CURL_ROOT=self.serve)
        r = subprocess.run(["bash", DOWNLOAD, self.queue, self.stage], capture_output=True,
                           text=True, env=env)
        self.assertEqual(r.returncode, 2)
        self.assertIn("header must name", r.stdout + r.stderr)

    def test_bad_queue_rows_stop_before_any_download(self):
        good = self.publish("GCA_000000001.1", ">c\nACGT\n")
        for row, frag in ((("../outside", "GCA_000000001.1", good), "species label"),
                          (("alpha", "GCA_1.1", good), "accession"),
                          (("alpha", "GCA_000000002.1", good), "does not name"),
                          (("alpha", "GCA_000000001.1", "file:///etc/passwd"), "ftp_path"),
                          (("alpha", "GCA_000000002.1", "https:///GCA_000000002.1"), "ftp_path")):
            with self.subTest(frag):
                r = self.run_dl([("alpha", "GCA_000000001.1", good), row]
                                if row[1] != "GCA_000000001.1" else [row])
                self.assertEqual(r.returncode, 2, r.stdout + r.stderr)
                self.assertIn(frag, r.stderr)
                self.assertFalse(os.path.exists(self.log))          # nothing fetched
                self.assertFalse(os.path.exists(os.path.join(self.tmp, "outside")))
        r = self.run_dl([("alpha", "GCA_000000001.1", good)] * 2)
        self.assertEqual(r.returncode, 2)
        self.assertIn("duplicate accession", r.stderr)

    def test_failed_final_move_leaves_neither_file(self):
        with open(os.path.join(self.bin, "mv"), "w") as fh:
            fh.write('#!/bin/sh\ncase "$2" in *.fna.ok) exit 1;; esac\nexec /bin/mv "$@"\n')
        os.chmod(os.path.join(self.bin, "mv"), 0o755)
        url = self.publish("GCA_000000001.1", ">c\nACGT\n")
        r = self.run_dl([("alpha", "GCA_000000001.1", url)])
        self.assertEqual(r.returncode, 1, r.stdout)
        d = os.path.join(self.stage, "alpha")
        self.assertEqual([x for x in os.listdir(d) if not x.startswith(".")], [])

    def test_gap_report_queue_downloads_as_is(self):
        # end to end: 250 gaps -> corpus_download_queue.tsv -> this tool
        t = Tree(os.path.join(self.tmp, "tree"))
        o1 = os.path.join(self.tmp, "b")
        r = run("build", "--sources", t.sources, "--out", o1)
        self.assertEqual(r.returncode, 0, r.stderr)
        corpus = os.path.join(o1, "corpus.tsv")
        summ = os.path.join(self.tmp, "summary.txt")
        missing = ["GCA_000000010.1", "GCA_000000011.1"]
        with open(summ, "w") as fh:
            fh.write(SUMMARY_HDR + "\n")
            for a in missing:
                fh.write(summary_line(a, 101) + "\n")
            fh.write(summary_line("GCA_000000002.3", 101) + "\n")    # version_mismatch
        import gzip
        for a in missing:                              # serve https://stub/<acc>/
            d = os.path.join(self.serve, "stub", a)
            os.makedirs(d)
            gz = os.path.join(d, a + "_genomic.fna.gz")
            with gzip.open(gz, "wt") as fh:
                fh.write(">%s\nACGT\n" % a)
            with open(os.path.join(d, "md5checksums.txt"), "w") as fh:
                fh.write("%s  ./%s_genomic.fna.gz\n" % (hashlib.md5(open(gz, "rb").read()).hexdigest(), a))
        o2 = os.path.join(self.tmp, "g")
        r = run("gaps", "--sources", t.sources, "--corpus", corpus, "--corpus-sha256", sha(corpus),
                "--summary", summ, "--summary-sha256", sha(summ), "--cutoff", "2026-09-01",
                "--out", o2)
        self.assertEqual(r.returncode, 0, r.stderr)
        env = dict(os.environ, PATH=self.bin + os.pathsep + os.environ["PATH"],
                   CURL_LOG=self.log, CURL_ROOT=self.serve, TRIES="1", BACKOFF="0")
        r = subprocess.run(["bash", DOWNLOAD, os.path.join(o2, "corpus_download_queue.tsv"),
                            self.stage], capture_output=True, text=True, env=env)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertIn("COMPLETE", r.stdout)
        self.assertEqual(sorted(os.listdir(os.path.join(self.stage, "alpha"))),
                         sorted([a + ".fna" for a in missing] + [a + ".fna.ok" for a in missing]))

if __name__ == "__main__":
    unittest.main(verbosity=2)
