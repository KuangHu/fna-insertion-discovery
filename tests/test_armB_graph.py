#!/usr/bin/env python3
"""End-to-end tests of scripts/30_armB_graph.py against SPEC item M1-B (B-1..B-10).

minigraph and gfatools are replaced by stub programs driven by a JSON spec, so
the tests are hermetic and run on the login node in seconds. Genomes are
synthetic and built so every expected coordinate is known by construction:

    contig = PADL + L + <allele> + R + PADR          (allele = '' for empty)

The allele interior of a '+' genome starts at p = len(PADL) + len(L); a '-'
genome stores the reverse complement of the whole contig.
"""
import fcntl
import hashlib
import json
import os
import random
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SCRIPT = os.path.join(ROOT, "scripts", "30_armB_graph.py")
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from lib.fasta import Fasta, revcomp      # noqa: E402

RNG = random.Random(11)


def rseq(n):
    return "".join(RNG.choice("ACGT") for _ in range(n))


PADL, L, R, PADR = rseq(600), rseq(700), rseq(700), rseq(600)
INS, INS2, INT, BIG = rseq(800), rseq(2000), rseq(300), rseq(6000)
INS_TSD = R[:9] + rseq(791)          # begins with the 9 bp that follow the site
X40, X80 = rseq(40), rseq(80)
SEGS = {"s1": L, "s2": INS, "s3": R, "s4": INT, "s5": INS2, "s6": BIG, "s7": INS_TSD,
        "s8": X40, "s9": X80}
P = len(PADL) + len(L)                # allele interior start on a '+' contig
WALK = {"": "*", INS: ">s2", INS2: ">s5", INT: ">s4", BIG: ">s6", INS_TSD: ">s7",
        X40: ">s8", X80: ">s9"}

MINIGRAPH_STUB = r'''
import json, os, sys
spec = json.load(open(os.environ["ARMB_STUB_SPEC"]))
if sys.argv[1:] == ["--version"]:
    print(spec.get("mg_version", "0.0-stub"))
    sys.exit(0)
with open(spec["log"], "a") as f:
    f.write(" ".join(sys.argv[1:]) + "\n")
if "-cxggs" in sys.argv:
    if spec.get("fail_build"):
        sys.exit(1)
    sys.stdout.write(spec["gfa"])
elif "--call" in sys.argv:
    fa = os.path.basename(sys.argv[-1])
    sys.stdout.write("".join(l + "\n" for l in spec["calls"].get(fa, [])))
else:
    sys.exit(2)
'''
SAMTOOLS_STUB = r'''
import os, sys
if sys.argv[1:] == ["--version"]:
    print("samtools 0.0-stub")
    sys.exit(0)
assert sys.argv[1] == "faidx"
path = sys.argv[2]
data = open(path, "rb").read()
out, i, n = [], 0, len(data)
while i < n:
    j = data.index(b"\n", i)
    name = data[i + 1:j].split()[0].decode()
    off = j + 1
    k = data.find(b"\n>", off)
    end = n if k == -1 else k + 1
    lines = data[off:end].split(b"\n")
    lines = [l for l in lines if l]
    ln = sum(len(l) for l in lines)
    lb = len(lines[0]) if lines else 0
    out.append("%s\t%d\t%d\t%d\t%d\n" % (name, ln, off, lb, lb + 1))
    i = end
open(path + ".fai", "w").writelines(out)
'''
GFATOOLS_STUB = r'''
import json, os, sys
spec = json.load(open(os.environ["ARMB_STUB_SPEC"]))
if sys.argv[1:] == ["version"]:
    print("gfatools: 0.0-stub")
    sys.exit(0)
if spec.get("fail_bubble"):
    sys.exit(1)
sys.stdout.write("".join(l + "\n" for l in spec["bubbles"]))
'''


def write_fasta(path, records, width=60, eol="\n"):
    """FASTA plus a correct .fai, so the tests never need samtools."""
    fai, off = [], 0
    with open(path, "w", newline="") as fh:
        for name, seq in records.items():
            head = ">" + name + eol
            fh.write(head)
            off += len(head)
            fai.append("%s\t%d\t%d\t%d\t%d\n" % (name, len(seq), off, width, width + len(eol)))
            for i in range(0, len(seq), width):
                line = seq[i:i + width] + eol
                fh.write(line)
                off += len(line)
    with open(path + ".fai", "w") as fh:
        fh.writelines(fai)


def bubble_line(key=("bk", 700, 1500)):
    return "\t".join(map(str, [key[0], key[1], key[2], 3, 2, 0, 0, 800, -1, -1, -1,
                               "s1,s2,s3", "*", INS]))


class Panel:
    """A temp dir with genomes, stub tools, a manifest and a run helper."""

    def __init__(self, root):
        self.root = root
        self.spec = {"gfa": "".join("S\t%s\t%s\n" % kv for kv in SEGS.items()),
                     "bubbles": [bubble_line()], "calls": {},
                     "log": os.path.join(root, "tool.log")}
        self.paths = []
        self.bin = os.path.join(root, "bin")
        os.makedirs(self.bin)
        for name, body in (("minigraph", MINIGRAPH_STUB), ("gfatools", GFATOOLS_STUB),
                           ("samtools", SAMTOOLS_STUB)):
            p = os.path.join(self.bin, name)
            with open(p, "w") as fh:
                fh.write("#!%s\n%s" % (sys.executable, body))
            os.chmod(p, 0o755)
        self.out = os.path.join(root, "out")

    def genome(self, name, allele, strand="+", snp_in_context=False, key=("bk", 700, 1500),
               call=None, subdir="g", contig="c1", extra_copy=False, pad_extra=0):
        tail = PADR + "A" * pad_extra
        if extra_copy:      # an exact copy of the empty-site needle, ~300 bp downstream
            tail = PADR[:300] + L[-50:] + allele + R[:50] + PADR[300:]
        seq = PADL + L + allele + R + tail
        if snp_in_context:                       # one SNP 10 bp left of the site
            i = P - 10
            seq = seq[:i] + ("A" if seq[i] != "A" else "C") + seq[i + 1:]
        if strand == "-":
            seq = revcomp(seq)
            T = len(seq)
            qs, qe = T - P - len(allele) - 5, T - P + 7       # deliberate slack
        else:
            qs, qe = P - 7, P + len(allele) + 5
        d = os.path.join(self.root, subdir)
        os.makedirs(d, exist_ok=True)
        path = os.path.join(d, name + ".fna")
        write_fasta(path, {contig: seq})
        self.paths.append(path)
        rec = call if call is not None else "%s:%d:%s:%s:%d:%d" % (
            WALK[allele], len(allele), strand, contig, qs, qe)
        self.spec["calls"].setdefault(name + ".fna", []).append(
            "%s\t%d\t%d\t>s1\t>s3\t%s" % (key[0], key[1], key[2], rec))
        return seq

    def run(self, *extra, env_commit="0" * 40):
        with open(os.path.join(self.root, "manifest"), "w") as fh:
            fh.write("".join(p + "\n" for p in self.paths))
        spec_path = os.path.join(self.root, "spec.json")
        with open(spec_path, "w") as fh:
            json.dump(self.spec, fh)
        env = dict(os.environ, ARMB_STUB_SPEC=spec_path)
        env.pop("GIT_COMMIT", None)
        if env_commit:
            env["GIT_COMMIT"] = env_commit
        cmd = [sys.executable, "-I", SCRIPT, "--manifest", os.path.join(self.root, "manifest"),
               "--outdir", self.out, "--minigraph", os.path.join(self.bin, "minigraph"),
               "--gfatools", os.path.join(self.bin, "gfatools"),
               "--samtools", os.path.join(self.bin, "samtools"), "--threads", "1",
               "--min-genome-bytes", "1", "--min-genome-frac", "0"] + list(extra)
        return subprocess.run(cmd, capture_output=True, text=True, env=env)

    def events(self, name="armB_events.tsv"):
        with open(os.path.join(self.out, name)) as fh:
            head = fh.readline().rstrip("\n").split("\t")
            return [dict(zip(head, l.rstrip("\n").split("\t"))) for l in fh]

    def marker(self):
        with open(os.path.join(self.out, "armB_complete.json")) as fh:
            return json.load(fh)

    def tool_calls(self):
        if not os.path.exists(self.spec["log"]):
            return []
        with open(self.spec["log"]) as fh:
            return fh.read().splitlines()


class TArmB(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="tarmb.", dir=os.environ.get("TMPDIR"))
        self.p = Panel(self.tmp)

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def ok(self, r):
        self.assertEqual(r.returncode, 0, r.stderr)

    def fails(self, r, code, frag):
        self.assertEqual(r.returncode, code, r.stderr)
        self.assertIn(frag, r.stderr)
        for name in ("armB_events.tsv", "armB_large_events.tsv", "armB_inserts.fna",
                     "armB_complete.json"):
            self.assertFalse(os.path.exists(os.path.join(self.p.out, name)), name)
        self.assertFalse(any(d.startswith(".armB_tmp") for d in os.listdir(self.p.out)))

    # ---------------- B-5 / B-6 positive cases
    def test_verified_plus_frame_and_published_outputs(self):
        self.p.genome("E", "")
        self.p.genome("F", INS)
        self.ok(self.p.run())
        ev = self.p.events()
        self.assertEqual(len(ev), 1)
        e = ev[0]
        self.assertEqual(e["frame_status"], "verified")
        self.assertEqual((int(e["junction_L"]), int(e["junction_R"]), int(e["junction_span"])),
                         (P, P, 0))
        self.assertEqual((e["empty_ref_genome"], e["filled_ref_genome"]), ("E", "F"))
        self.assertEqual(int(e["insert_len_direct"]), 800)
        self.assertEqual(e["insert_md5"], hashlib.md5(INS.encode()).hexdigest())
        self.assertEqual((e["n_empty"], e["n_filled"], e["n_intermediate"], e["n_missing"]),
                         ("1", "1", "0", "0"))
        self.assertEqual(e["flags"], ".")
        m = self.p.marker()
        self.assertEqual(m["status"], "complete")
        self.assertEqual(m["git_commit_submit"], "0" * 40)
        for name, sha in m["outputs_sha256"].items():
            with open(os.path.join(self.p.out, name), "rb") as fh:
                self.assertEqual(hashlib.sha256(fh.read()).hexdigest(), sha, name)
        self.assertFalse(any(d.startswith(".armB_tmp") for d in os.listdir(self.p.out)))
        with open(os.path.join(self.p.out, "armB_empty_sites.fna")) as fh:
            site = fh.read().split("\n")[1]
        self.assertEqual(site, L[-500:] + R[:500])

    def test_minus_strand_frame(self):
        seq = self.p.genome("E", "", strand="-")
        self.p.genome("F", INS)
        self.ok(self.p.run())
        e = self.p.events()[0]
        T = len(seq)
        self.assertEqual(e["frame_status"], "verified")
        self.assertEqual(e["empty_ref_strand"], "-")
        self.assertEqual((int(e["junction_L"]), int(e["junction_R"])), (T - P, T - P))
        self.assertEqual(e["insert_md5"], hashlib.md5(INS.encode()).hexdigest())

    def test_minus_strand_with_target_span(self):
        # empty allele keeps 300 bp that the filled allele replaces: span 300
        seq = self.p.genome("E", INT, strand="-")
        self.p.genome("F", INS2)
        self.ok(self.p.run())
        e = self.p.events()[0]
        T = len(seq)
        self.assertEqual(e["frame_status"], "verified")
        self.assertEqual((int(e["junction_L"]), int(e["junction_R"])), (T - P - 300, T - P))
        self.assertEqual(int(e["junction_span"]), 300)
        self.assertEqual(e["architecture"], "asymmetric_candidate")

    def test_frame_mismatch_gives_no_coordinates(self):
        self.p.genome("E", "", snp_in_context=True)
        self.p.genome("F", INS)
        self.ok(self.p.run())
        e = self.p.events()[0]
        self.assertEqual(e["frame_status"], "unverified")
        self.assertEqual((e["junction_L"], e["junction_R"]), ("-1", "-1"))
        self.assertEqual(e["architecture"], "no_empty_frame")
        self.assertIn("frame_unverified", e["flags"])
        self.assertEqual((int(e["empty_ref_call_start"]), int(e["empty_ref_call_end"])),
                         (P - 7, P + 5))
        self.assertEqual(e["insert_md5"], hashlib.md5(INS.encode()).hexdigest())

    def test_neighbouring_copy_is_not_the_frame(self):
        # true site carries a SNP; an exact copy lies outside the call interval
        self.p.genome("E", "", snp_in_context=True, extra_copy=True)
        self.p.genome("F", INS)
        self.ok(self.p.run())
        e = self.p.events()[0]
        self.assertEqual(e["frame_status"], "off_call_interval")
        self.assertEqual((e["junction_L"], e["junction_R"]), ("-1", "-1"))

    def test_copy_inside_a_wide_call_interval_is_not_the_frame(self):
        # true site has a SNP; the exact copy (~1,000 bp downstream) lies INSIDE
        # a wide call interval, so overlap alone would accept it
        wide = "*:0:+:c1:%d:%d" % (P - 7, P + 1150)
        self.p.genome("E", "", snp_in_context=True, extra_copy=True, call=wide)
        self.p.genome("F", INS)
        self.ok(self.p.run())
        e = self.p.events()[0]
        self.assertEqual(e["frame_status"], "interval_mismatch")
        self.assertEqual((e["junction_L"], e["junction_R"]), ("-1", "-1"))

    def test_two_exact_copies_are_ambiguous(self):
        self.p.genome("E", "", extra_copy=True)
        self.p.genome("F", INS)
        self.ok(self.p.run())
        self.assertEqual(self.p.events()[0]["frame_status"], "ambiguous")

    def test_frame_reason_belongs_to_the_declared_reference(self):
        # E1 (most contiguous: longest) has no exact copy; E2 has two copies
        self.p.genome("E1", "", snp_in_context=True, pad_extra=500)
        self.p.genome("E2", "", extra_copy=True)
        self.p.genome("F", INS)
        self.ok(self.p.run())
        e = self.p.events()[0]
        self.assertEqual((e["empty_ref_genome"], e["frame_status"]), ("E1", "unverified"))

    def test_prefers_an_empty_genome_that_verifies(self):
        self.p.genome("E1", "", snp_in_context=True)
        self.p.genome("E2", "")
        self.p.genome("F", INS)
        self.ok(self.p.run())
        e = self.p.events()[0]
        self.assertEqual((e["frame_status"], e["empty_ref_genome"]), ("verified", "E2"))

    # ---------------- B-3 / B-5 classification and consistency
    def test_intermediate_and_missing_counted(self):
        self.p.genome("E", "")
        self.p.genome("F", INS)
        self.p.genome("I", INT)
        self.p.genome("M", "", call="*")
        self.ok(self.p.run())
        e = self.p.events()[0]
        self.assertEqual((e["n_empty"], e["n_filled"], e["n_intermediate"], e["n_missing"]),
                         ("1", "1", "1", "1"))
        self.assertIn("intermediate_1", e["flags"])
        self.assertIn("missing_in_1", e["flags"])

    def test_multiallelic_and_length_mismatch_flagged(self):
        self.p.genome("E", "")
        self.p.genome("F1", INS)
        self.p.genome("F2", INS2)
        self.ok(self.p.run())
        e = self.p.events()[0]
        self.assertEqual(int(e["insert_size"]), 1400)          # median(800, 2000)
        self.assertEqual(e["n_length_classes"], "3")
        self.assertIn("multiallelic", e["flags"])
        self.assertIn("insert_len_mismatch", e["flags"])
        self.assertIn(int(e["insert_len_direct"]), (800, 2000))

    def test_intermediate_lengths_cannot_chain_classes(self):
        self.p.genome("E", "")
        self.p.genome("E2", X40)          # empty (within 50 of the shortest)
        self.p.genome("I", X80)           # intermediate
        self.p.genome("F", INS)
        self.ok(self.p.run())
        e = self.p.events()[0]
        self.assertEqual((e["n_empty"], e["n_intermediate"], e["n_filled"]), ("2", "1", "1"))
        self.assertEqual(e["n_length_classes"], "3")
        self.assertIn("multiallelic", e["flags"])

    def test_g4_provenance_in_marker(self):
        self.p.genome("E", "")
        self.p.genome("F", INS)
        self.ok(self.p.run())
        m = self.p.marker()
        self.assertEqual(m["tool_versions"]["minigraph"], "0.0-stub")
        self.assertEqual(m["tool_versions"]["gfatools"], "gfatools: 0.0-stub")
        for k in ("host", "slurm_job_id", "slurm_array_task_id", "git_commit_runtime",
                  "git_dirty", "params", "argv"):
            self.assertIn(k, m)
        for rec, path in zip(m["inputs"], self.p.paths):
            with open(path, "rb") as fh:
                self.assertEqual(rec["sha256"], hashlib.sha256(fh.read()).hexdigest())

    def test_tsd_is_annotation_only(self):
        self.p.genome("E", "")
        self.p.genome("F", INS_TSD)
        self.ok(self.p.run())
        e = self.p.events()[0]
        self.assertEqual((e["tsd_len"], e["tsd_seq"], e["tsd_side"]), ("9", R[:9], "right"))
        # the event without a TSD (previous tests) is emitted just the same
        self.assertEqual(len(self.p.events()), 1)

    # ---------------- B-4 size gate
    def test_large_event_parked_with_sequence(self):
        self.p.genome("E", "")
        self.p.genome("F", BIG)
        self.ok(self.p.run())
        self.assertEqual(self.p.events(), [])
        lg = self.p.events("armB_large_events.tsv")
        self.assertEqual(len(lg), 1)
        self.assertEqual(int(lg[0]["insert_len_direct"]), 6000)
        with open(os.path.join(self.p.out, "armB_large_inserts.fna")) as fh:
            self.assertEqual(fh.read().split("\n")[1], BIG)
        self.assertEqual(self.p.marker()["counts"]["parked"], 1)
        self.assertNotIn("parent_large_event_id", open(os.path.join(
            self.p.out, "armB_events.tsv")).readline())

    def test_below_min_insert_is_counted(self):
        self.p.genome("E", "")
        self.p.genome("S", INT)
        self.ok(self.p.run())
        self.assertEqual(self.p.events(), [])
        self.assertEqual(self.p.marker()["counts"]["below_min_insert"], 1)

    # ---------------- B-1 failure handling
    def test_tool_failure_leaves_no_output(self):
        self.p.genome("E", "")
        self.p.genome("F", INS)
        os.makedirs(self.p.out)
        with open(os.path.join(self.p.out, "armB_events.tsv"), "w") as fh:
            fh.write("stale\tfile\nfrom\tbefore\n")
        self.p.spec["fail_build"] = True
        self.fails(self.p.run(), 2, "external tool failed")

    def test_call_contig_not_in_fasta_is_fatal(self):
        self.p.genome("E", "")
        self.p.genome("F", INS, call=">s2:800:+:cX:593:1405")   # contig not in FASTA
        self.fails(self.p.run(), 2, "is not in F's FASTA")

    def test_crash_midway_leaves_no_final_file(self):
        """An exception while events are being written (here injected into the
        TSD annotation, after the temporary outputs are open and partly
        written) must leave no final output and no marker."""
        import importlib.util
        self.p.genome("E", "")
        self.p.genome("F", INS)
        self.p.run("--allow-zero-bubbles")          # writes manifest/spec; a good run
        self.assertTrue(os.path.exists(os.path.join(self.p.out, "armB_complete.json")))
        spec = importlib.util.spec_from_file_location("armb_crash", SCRIPT)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)

        def boom(*a, **k):
            raise RuntimeError("injected crash mid-write")
        mod.find_tsd = boom
        old = dict(os.environ)
        os.environ.update(ARMB_STUB_SPEC=os.path.join(self.tmp, "spec.json"), GIT_COMMIT="0" * 40)
        try:
            rc = mod.main(["--manifest", os.path.join(self.tmp, "manifest"), "--outdir", self.p.out,
                           "--minigraph", os.path.join(self.p.bin, "minigraph"),
                           "--gfatools", os.path.join(self.p.bin, "gfatools"), "--threads", "1",
                           "--min-genome-bytes", "1", "--min-genome-frac", "0"])
        finally:
            os.environ.clear()
            os.environ.update(old)
        self.assertEqual(rc, 2)
        for name in ("armB_events.tsv", "armB_inserts.fna", "armB_complete.json"):
            self.assertFalse(os.path.exists(os.path.join(self.p.out, name)), name)
        self.assertFalse(any(d.startswith(".armB_tmp") for d in os.listdir(self.p.out)))

    def test_unreadable_fasta_is_fatal(self):
        self.p.genome("E", "")
        self.p.genome("F", INS)
        with open(self.p.paths[1] + ".fai", "w") as fh:
            fh.write("c1\tnot_a_number\n")
        self.fails(self.p.run(), 2, "cannot index/read")

    def test_reuse_is_bound_to_graph_and_inputs(self):
        self.p.genome("E", "")
        self.p.genome("F", INS)
        self.ok(self.p.run("--reuse"))
        n1 = len(self.p.tool_calls())                  # 1 build + 2 calls
        self.assertEqual(n1, 3)
        self.ok(self.p.run("--reuse"))
        self.assertEqual(len(self.p.tool_calls()), n1)   # everything reused
        # graph rebuilt with different content: old call files must not be reused
        os.remove(os.path.join(self.p.out, "graph.gfa"))
        self.p.spec["gfa"] += "S\ts99\tACGT\n"
        self.ok(self.p.run("--reuse"))
        new = self.p.tool_calls()[n1:]
        self.assertEqual(sum("-cxggs" in c for c in new), 1)
        self.assertEqual(sum("--call" in c for c in new), 2)

    def test_reuse_ignores_an_mtime_only_change(self):
        # identity is content (sha256): touching a file does not invalidate it
        self.p.genome("E", "")
        self.p.genome("F", INS)
        self.ok(self.p.run("--reuse"))
        n1 = len(self.p.tool_calls())
        st = os.stat(self.p.paths[0])
        os.utime(self.p.paths[0], ns=(st.st_atime_ns, st.st_mtime_ns + 10**9))
        self.ok(self.p.run("--reuse"))
        self.assertEqual(self.p.tool_calls()[n1:], [])

    def test_stale_marker_removed_before_outputs(self):
        self.p.genome("E", "")
        self.p.genome("F", INS)
        self.ok(self.p.run())
        self.p.spec["fail_build"] = True
        self.fails(self.p.run(), 2, "external tool failed")   # marker + outputs gone

    def test_git_commit_required(self):
        self.p.genome("E", "")
        self.p.genome("F", INS)
        self.fails(self.p.run(env_commit=None), 2, "GIT_COMMIT")

    def test_missing_executable_reports_fatal(self):
        self.p.genome("E", "")
        self.p.genome("F", INS)
        self.fails(self.p.run("--minigraph", "/nonexistent/minigraph"), 2, "FATAL: FileNotFoundError")

    def test_stale_fai_is_rejected(self):
        self.p.genome("E", "")
        self.p.genome("F", INS)
        path = self.p.paths[1]
        with open(path) as fh:
            data = fh.read()
        with open(path, "w") as fh:          # header grows by exactly one sequence line
            fh.write(data.replace(">c1\n", ">c1 " + "x" * 60 + "\n", 1))
        self.fails(self.p.run(), 2, "stale or wrong .fai")

    def test_call_file_must_cover_every_bubble(self):
        self.p.genome("E", "")
        self.p.genome("F", INS)
        self.p.spec["calls"]["F.fna"] = []                    # tool exits 0, empty file
        self.fails(self.p.run(), 2, "does not cover the graph's bubbles")

    def test_concurrent_run_refused(self):
        self.p.genome("E", "")
        self.p.genome("F", INS)
        os.makedirs(self.p.out)
        with open(os.path.join(self.p.out, ".armB.lock"), "w") as lk:
            fcntl.flock(lk, fcntl.LOCK_EX | fcntl.LOCK_NB)
            r = self.p.run()
        self.assertEqual(r.returncode, 2, r.stderr)
        self.assertIn("another Arm B run", r.stderr)
        self.ok(self.p.run())                       # lock released: runs normally

    def test_refused_run_writes_no_failure_record(self):
        self.p.genome("E", "")
        self.p.genome("F", INS)
        self.ok(self.p.run())
        with open(os.path.join(self.p.out, ".armB.lock"), "w") as lk:
            fcntl.flock(lk, fcntl.LOCK_EX | fcntl.LOCK_NB)
            r = self.p.run()
        self.assertEqual(r.returncode, 2)
        self.assertFalse(os.path.exists(os.path.join(self.p.out, "armB_failed.json")))
        self.assertTrue(os.path.exists(os.path.join(self.p.out, "armB_complete.json")))

    def test_reuse_bound_to_content_not_mtime(self):
        self.p.genome("E", "")
        self.p.genome("F", INS)
        self.ok(self.p.run("--reuse"))
        n1 = len(self.p.tool_calls())
        path = self.p.paths[1]
        st = os.stat(path)
        with open(path, "r+") as fh:                  # same size, same mtime, new content
            fh.seek(st.st_size - 10)
            c = fh.read(1)
            fh.seek(st.st_size - 10)
            fh.write("A" if c != "A" else "C")
        os.utime(path, ns=(st.st_atime_ns, st.st_mtime_ns))
        self.ok(self.p.run("--reuse"))
        new = self.p.tool_calls()[n1:]
        # the graph is rebuilt (an input changed); F's calls are redone. The stub
        # graph comes out byte-identical, so E's call file (same graph sha256,
        # same E content) is legitimately reused.
        self.assertEqual(sum("-cxggs" in c for c in new), 1)
        self.assertEqual([c.split()[-1].endswith("F.fna") for c in new if "--call" in c], [True])

    def test_reuse_bound_to_generator_version_executable_and_args(self):
        self.p.genome("E", "")
        self.p.genome("F", INS)
        self.ok(self.p.run("--reuse"))
        first = self.p.marker()["cache"]
        self.assertFalse(first["graph"]["reused"])
        self.assertEqual(first["graph"]["meta"]["generator"]["version"], "0.0-stub")
        self.assertEqual(first["graph"]["meta"]["generator"]["args"], ["-cxggs", "-t", "1"])

        def redone(n):
            new = self.p.tool_calls()[n:]
            return sum("-cxggs" in c for c in new), sum("--call" in c for c in new)
        # another minigraph version; the stub graph comes out byte-identical,
        # so the call files are invalidated by the generator, not by the graph
        n = len(self.p.tool_calls())
        self.p.spec["mg_version"] = "0.1-stub"
        self.ok(self.p.run("--reuse"))
        self.assertEqual(redone(n), (1, 2))
        m = self.p.marker()
        self.assertEqual(m["tool_versions"]["minigraph"], "0.1-stub")
        self.assertEqual(m["cache"]["graph"]["meta"]["generator"]["version"], "0.1-stub")
        self.assertEqual({c["meta"]["generator"]["version"] for c in m["cache"]["calls"].values()},
                         {"0.1-stub"})
        # same version string, different executable
        n = len(self.p.tool_calls())
        with open(os.path.join(self.p.bin, "minigraph"), "a") as fh:
            fh.write("\n# rebuilt\n")
        self.ok(self.p.run("--reuse"))
        self.assertEqual(redone(n), (1, 2))
        # different arguments
        n = len(self.p.tool_calls())
        self.ok(self.p.run("--reuse", "--threads", "2"))
        self.assertEqual(redone(n), (1, 2))
        # nothing changed: everything reused
        n = len(self.p.tool_calls())
        self.ok(self.p.run("--reuse", "--threads", "2"))
        self.assertEqual(redone(n), (0, 0))

    def test_cache_without_generator_record_is_not_reused(self):
        # metadata written before generators were recorded cannot be trusted
        self.p.genome("E", "")
        self.p.genome("F", INS)
        self.ok(self.p.run("--reuse"))
        n = len(self.p.tool_calls())
        for rel in ("graph.gfa.meta.json", "calls/E.bed.meta.json"):
            path = os.path.join(self.p.out, rel)
            with open(path) as fh:
                meta = json.load(fh)
            del meta["generator"]
            with open(path, "w") as fh:
                json.dump(meta, fh)
        self.ok(self.p.run("--reuse"))
        new = self.p.tool_calls()[n:]
        self.assertEqual(sum("-cxggs" in c for c in new), 1)
        # the rebuilt graph is byte-identical and F's meta is intact: only E is redone
        self.assertEqual([c.split()[-1].endswith("E.fna") for c in new if "--call" in c], [True])
        calls = self.p.marker()["cache"]["calls"]
        self.assertEqual((calls["E"]["reused"], calls["F"]["reused"]), (False, True))

    def test_marker_carries_original_record_of_reused_files(self):
        self.p.genome("E", "")
        self.p.genome("F", INS)
        self.ok(self.p.run("--reuse", env_commit="a" * 40))
        self.ok(self.p.run("--reuse", env_commit="b" * 40))
        m = self.p.marker()
        self.assertEqual(m["git_commit_submit"], "b" * 40)
        self.assertTrue(m["cache"]["graph"]["reused"])
        self.assertEqual(m["cache"]["graph"]["meta"]["created"]["git_commit_submit"], "a" * 40)
        for sn in ("E", "F"):
            c = m["cache"]["calls"][sn]
            self.assertTrue(c["reused"])
            self.assertEqual(c["meta"]["created"]["git_commit_submit"], "a" * 40)
            self.assertEqual(c["meta"]["generator"]["args"], ["-cxasm", "--call", "-t", "1"])

    def test_marker_invalidated_before_startup_checks(self):
        self.p.genome("E", "")
        self.p.genome("F", INS)
        self.ok(self.p.run())
        self.fails(self.p.run("--len-tol", "600"), 2, "--len-tol < --min-insert")
        self.ok(self.p.run())
        self.fails(self.p.run(env_commit=None), 2, "GIT_COMMIT")

    def test_call_past_contig_end_is_fatal(self):
        self.p.genome("E", "")
        self.p.genome("F", INS, call=">s2:800:+:c1:593:999999")
        self.fails(self.p.run(), 2, "runs past the end")

    def test_unknown_node_in_unused_intermediate_walk_is_fatal(self):
        self.p.genome("E", "")
        self.p.genome("F", INS)
        self.p.genome("I", INT, call=">s404:300:+:c1:1293:1598")
        self.fails(self.p.run(), 2, "unknown graph node")

    def test_walk_length_must_match_call_for_every_genome(self):
        self.p.genome("E", "")
        self.p.genome("F", INS)
        self.p.genome("I", INT, call=">s4:310:+:c1:1293:1598")      # INT is 300 bp
        self.fails(self.p.run(), 2, "has 300 bp but the call says 310")

    # ---------------- B-2 parsing
    def test_bubble_columns_8_to_13_checked(self):
        for bad in (bubble_line().replace("\t800\t", "\tx800\t", 1),
                    bubble_line().replace(INS, INS[:-1]),
                    bubble_line().replace("\ts1,s2,s3\t", "\ts1,,s3\t"),
                    bubble_line().replace("\t*\t", "\tA\t")):
            self.p.spec["bubbles"] = [bad]
            if not self.p.paths:
                self.p.genome("E", "")
                self.p.genome("F", INS)
            self.fails(self.p.run(), 2, "unparseable gfatools bubble")

    def test_failure_summary_records_parse_counts(self):
        self.p.genome("E", "")
        self.p.genome("F", INS, call="weird:record")
        self.fails(self.p.run(), 2, "unparseable minigraph --call")
        with open(os.path.join(self.p.out, "armB_failed.json")) as fh:
            rec = json.load(fh)
        self.assertEqual((rec["status"], rec["details"]["unparseable_call_lines"]), ("failed", 1))
        self.ok(self._fix_and_rerun())
        self.assertFalse(os.path.exists(os.path.join(self.p.out, "armB_failed.json")))
        c = self.p.marker()["counts"]
        self.assertEqual((c["unparseable_call_lines"], c["unparseable_bubble_lines"]), (0, 0))

    def _fix_and_rerun(self):
        good = "bk\t700\t1500\t>s1\t>s3\t>s2:800:+:c1:%d:%d" % (P - 7, P + 805)
        self.p.spec["calls"]["F.fna"] = [good]
        return self.p.run()

    def test_samtools_version_recorded_when_it_builds_an_index(self):
        self.p.genome("E", "")
        self.p.genome("F", INS)
        self.ok(self.p.run())
        self.assertTrue(self.p.marker()["tool_versions"]["samtools"].startswith("not invoked"))
        os.remove(self.p.paths[1] + ".fai")
        self.ok(self.p.run())
        self.assertEqual(self.p.marker()["tool_versions"]["samtools"], "samtools 0.0-stub")

    def test_empty_sixth_column_is_fatal(self):
        self.p.genome("E", "")
        self.p.genome("F", INS, call="")
        self.fails(self.p.run(), 2, "unparseable minigraph --call")

    def test_duplicate_record_with_star_is_fatal(self):
        self.p.genome("E", "")
        self.p.genome("F", INS)
        self.p.spec["calls"]["F.fna"].append("bk\t700\t1500\t>s1\t>s3\t*")
        self.fails(self.p.run(), 2, "unparseable minigraph --call")

    def test_malformed_walk_is_fatal_even_when_not_used(self):
        self.p.genome("E", "")
        self.p.genome("F", INS)
        self.p.genome("I", INT, call="s4:300:+:c1:1293:1598")   # walk without orientation
        self.fails(self.p.run(), 2, "unparseable minigraph --call")

    def test_unparseable_bubble_line_is_fatal(self):
        self.p.genome("E", "")
        self.p.genome("F", INS)
        self.p.spec["bubbles"] = ["bk\t700\t1500\t3\t2\t0\t0\t800\t-1\t-1"]
        self.fails(self.p.run(), 2, "unparseable gfatools bubble")

    def test_unparseable_call_record_is_fatal(self):
        self.p.genome("E", "")
        self.p.genome("F", INS, call="weird:record")
        self.fails(self.p.run(), 2, "unparseable minigraph --call")

    # ---------------- inputs, parameters, zero bubbles
    def test_duplicate_sample_names_rejected(self):
        self.p.genome("E", "", subdir="a")
        self.p.genome("E", INS, subdir="b")
        self.fails(self.p.run(), 2, "duplicate names")

    def test_len_tol_must_be_below_min_insert(self):
        self.p.genome("E", "")
        self.p.genome("F", INS)
        self.fails(self.p.run("--len-tol", "600"), 2, "--len-tol < --min-insert")

    def test_zero_bubbles_fails_by_default(self):
        self.p.genome("E", "")
        self.p.genome("F", INS)
        self.p.spec["bubbles"] = []
        self.p.spec["calls"] = {k: [] for k in self.p.spec["calls"]}   # no bubbles, no records
        self.fails(self.p.run(), 3, "0 bubbles")
        self.ok(self.p.run("--allow-zero-bubbles"))
        self.assertEqual(self.p.marker()["counts"]["bubbles"], 0)


class TFasta(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="tfa.", dir=os.environ.get("TMPDIR"))

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def test_crlf_multiline_fetch(self):
        seq = rseq(1000)
        p = os.path.join(self.tmp, "x.fna")
        write_fasta(p, {"c": seq}, width=60, eol="\r\n")
        fa = Fasta(p)
        for a, b in ((0, 1000), (55, 185), (59, 61), (999, 1000), (120, 120)):
            self.assertEqual(fa.fetch("c", a, b), seq[a:b])

    def test_short_read_raises(self):
        seq = rseq(1000)
        p = os.path.join(self.tmp, "x.fna")
        write_fasta(p, {"c": seq})
        with open(p, "r+") as fh:
            fh.truncate(500)
        # caught when the index is bound to the file, or at the latest on fetch
        with self.assertRaises((IOError, ValueError)):
            Fasta(p).fetch("c", 0, 1000)

    def test_stale_index_rejected(self):
        seq = rseq(300)
        p = os.path.join(self.tmp, "x.fna")
        write_fasta(p, {"c": seq})
        with open(p) as fh:
            data = fh.read()
        with open(p, "w") as fh:
            fh.write(data.replace(">c\n", ">c " + "y" * 60 + "\n"))
        with self.assertRaises(ValueError):
            Fasta(p)

    def test_index_shorter_than_sequence_rejected(self):
        p = os.path.join(self.tmp, "x.fna")
        write_fasta(p, {"c": rseq(120)}, width=60)
        with open(p + ".fai", "w") as fh:
            fh.write("c\t60\t3\t60\t61\n")          # stale: declares 60 of 120 bp
        with self.assertRaises(ValueError):
            Fasta(p)

    def test_zero_length_index_entry_for_real_sequence_rejected(self):
        p = os.path.join(self.tmp, "x.fna")
        write_fasta(p, {"c": rseq(120)}, width=60)
        with open(p + ".fai", "w") as fh:
            fh.write("c\t0\t3\t60\t61\n")
        with self.assertRaises(ValueError):
            Fasta(p)

    def test_empty_or_truncated_index_rejected(self):
        p = os.path.join(self.tmp, "x.fna")
        write_fasta(p, {"a": rseq(130), "b": rseq(61)}, width=60)
        with open(p + ".fai") as fh:
            lines = fh.readlines()
        for content in ("", lines[0]):           # empty; second record missing
            with open(p + ".fai", "w") as fh:
                fh.write(content)
            with self.assertRaises(ValueError):
                Fasta(p)

    def test_iupac_reverse_complement(self):
        self.assertEqual(revcomp("AR"), "YT")
        self.assertEqual(revcomp("ACGTRYKMBVDHSWN"), "NWSDHBVKMRYACGT")
        self.assertEqual(revcomp("acgtrykmbvdhswn"), "nwsdhbvkmryacgt")

    def test_two_records_index_ok(self):
        p = os.path.join(self.tmp, "x.fna")
        a, b = rseq(130), rseq(61)
        write_fasta(p, {"a": a, "b": b}, width=60)
        fa = Fasta(p)
        self.assertEqual((fa.fetch("a", 50, 130), fa.fetch("b", 0, 61)), (a[50:], b))

    def raw(self, data, fai):
        p = os.path.join(self.tmp, "raw.fna")
        with open(p, "wb") as fh:
            fh.write(data)
        with open(p + ".fai", "w") as fh:
            fh.write(fai)
        return p

    def test_every_line_of_a_record_is_checked(self):
        # each passes the old first-line / last-base checks, and fetch() would
        # return a base of the right length from the wrong place
        cases = [
            (b">c\nACGT\nAC\nGTACGT\nACGT\n", "c\t16\t3\t4\t5\n"),      # review r07
            (b">c\nACGT\nACGTA\nCGT\nACGT\n", "c\t16\t3\t4\t5\n"),      # long, short
            (b">c\nACGT\r\nACGTC\nACGT\r\nAC\r\n", "c\t14\t3\t4\t6\n"),  # bare LF in CRLF
            (b">c\nACGT\nA>GT\nACGT\nAC\n", "c\t14\t3\t4\t5\n"),         # '>' inside
        ]
        for data, fai in cases:
            with self.subTest(data=data):
                with self.assertRaises(ValueError):
                    Fasta(self.raw(data, fai))

    def test_line_check_spans_read_chunks(self):
        # > 4 MB, so the record is checked in more than one read; the bad
        # line is far into the second chunk
        tab = bytes.maketrans(bytes(range(256)), b"ACGT" * 64)
        seq = os.urandom(4_500_000).translate(tab).decode()
        p = os.path.join(self.tmp, "big.fna")
        write_fasta(p, {"c": seq}, width=60)
        fa = Fasta(p)
        self.assertEqual(fa.fetch("c", 4_400_000, 4_400_130), seq[4_400_000:4_400_130])
        fa.close()
        with open(p, "r+b") as fh:          # move one line end by one byte, same size
            pos = 3 + 72_000 * 61 + 60
            fh.seek(pos)
            fh.write(b"A\n")
        with self.assertRaises(ValueError):
            Fasta(p)

    def test_every_open_rechecks_the_file(self):
        # no cache across opens: an in-place, same-size rewrite with the
        # timestamps restored is caught on the next open
        p = os.path.join(self.tmp, "x.fna")
        write_fasta(p, {"c": rseq(300)}, width=60)
        Fasta(p).close()
        st = os.stat(p)
        with open(p, "r+b") as fh:          # same size, mtime restored below
            fh.seek(3 + 60)
            fh.write(b"A\n")
        os.utime(p, ns=(st.st_atime_ns, st.st_mtime_ns))
        with self.assertRaises(ValueError):
            Fasta(p)

    def test_malformed_index_rejected(self):
        p = os.path.join(self.tmp, "x.fna")
        write_fasta(p, {"c": rseq(100)})
        with open(p + ".fai", "w") as fh:
            fh.write("c\t100\t3\t0\t1\n")              # line_bases 0
        with self.assertRaises(ValueError):
            Fasta(p)


if __name__ == "__main__":
    unittest.main(verbosity=2)
