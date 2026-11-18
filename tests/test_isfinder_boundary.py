#!/usr/bin/env python3
"""Tests of scripts/260_isfinder_boundary.py: the end projection, the classes,
and one end-to-end run on a made-up element (needs minimap2 on PATH)."""
import csv, importlib.util, json, os, random, shutil, subprocess, sys, tempfile, unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(HERE, "..", "scripts", "260_isfinder_boundary.py")
spec = importlib.util.spec_from_file_location("m260", SCRIPT)
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)


def rc(s):
    return s[::-1].translate(str.maketrans("ACGT", "TGCA"))


class TPlace(unittest.TestCase):
    def test_full_length_forward(self):
        # insert = 4 extra + IS(100) + 5 extra
        self.assertEqual(m.place(109, 4, 104, "+", 100, 0, 100), (4, 5))

    def test_clipped_ends_are_projected(self):
        # aligner dropped 3 IS bases at the start and 2 at the end
        self.assertEqual(m.place(109, 7, 102, "+", 100, 3, 98), (4, 5))

    def test_reverse_strand_swaps_the_tails(self):
        # IS reversed in the insert: the IS start lies at the insert's right
        self.assertEqual(m.place(109, 6, 101, "-", 100, 3, 98), (4, 5))

    def test_insert_ending_inside_the_element_is_negative(self):
        self.assertEqual(m.place(80, 0, 80, "+", 100, 20, 100), (-20, 0))

    def test_classes(self):
        self.assertEqual(m.classify(9, 0, 9, 3), ("exact", 0))
        self.assertEqual(m.classify(4, 5, 9, 3), ("exact", 0))
        self.assertEqual(m.classify(-1, 10, 9, 3), ("near", 0))   # one end inside the IS
        # repeat 9, overlap 10 by a chance base: cut rotated by one, or not
        self.assertEqual(m.classify(-1, 10, 10, 3), ("exact_shifted", -1))
        self.assertEqual(m.classify(0, 9, 10, 3), ("exact_shifted", -1))
        self.assertEqual(m.classify(-2, 11, 10, 3), ("near", -1))  # rotated past the slack
        self.assertEqual(m.classify(5, 5, 9, 3), ("near", 1))
        self.assertEqual(m.classify(300, 9, 9, 3), ("over", 300))
        self.assertEqual(m.classify(-200, 0, 9, 3), ("under", -209))


@unittest.skipUnless(shutil.which("minimap2"), "minimap2 not on PATH")
class TEndToEnd(unittest.TestCase):
    def test_run(self):
        rnd = random.Random(7)
        seq = lambda n: "".join(rnd.choice("ACGT") for _ in range(n))
        is1, is2 = seq(1200), seq(900)
        d = tempfile.mkdtemp(dir=os.environ.get("TMPDIR"))
        self.addCleanup(shutil.rmtree, d)
        isf = os.path.join(d, "IS.fna")
        with open(isf, "w") as fh:
            fh.write(">ISX_ISfam_grp\n%s\n>ISY_ISfam_grp\n%s\n" % (is1, is2))
        with open(os.path.join(d, "panel.tsv"), "w") as fh:
            fh.write("is_name\ttsd_literature\tnote\nISX\t9\t\n")
        rows = [  # id, insert, overlap, expected class
            ("E1", is1 + seq(9), 9, "exact"),
            ("E2", seq(4) + rc(is1) + seq(5), 9, "exact"),
            ("E3", is1 + seq(9), 0, "over"),
            ("E4", seq(400) + is1 + seq(9), 9, "over"),
            ("E5", is1[60:] + seq(9), 9, "under"),
            ("E6", seq(1500), 5, None),
            ("E7", is2 + seq(3), 3, "exact"),
        ]
        cols = ["event_id", "species", "insert_seq", "inserted_len", "junction_overlap_bp",
                "insert_frame_status", "decomposition_method", "S1_structurally_clean"]
        ev = os.path.join(d, "x_events.tsv")
        with open(ev, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n"); w.writerow(cols)
            for i, s, ov, _ in rows:
                w.writerow([i, "sp", s, len(s), ov, "coord", "exact", "1"])
        out = os.path.join(d, "out")
        cmd = [sys.executable, SCRIPT, "--isfinder", isf, "--isfinder-sha256", m.sha256_file(isf),
               "--events", ev, "--panel", os.path.join(d, "panel.tsv"), "--out", out,
               "--threads", "1", "--min-n", "1"]
        p = subprocess.run(cmd, capture_output=True, text=True)
        self.assertEqual(p.returncode, 0, p.stderr)
        got = {r["event_id"]: r for r in csv.DictReader(open(os.path.join(out, "per_event.tsv")), delimiter="\t")}
        for i, _, _, want in rows:
            if want is None:
                self.assertNotIn(i, got)
            else:
                self.assertEqual(got[i]["class"], want, i)
        self.assertEqual(got["E2"]["strand"], "-")
        self.assertEqual((got["E2"]["left_extra"], got["E2"]["right_extra"]), ("4", "5"))
        res = json.load(open(os.path.join(out, "result.json")))
        self.assertEqual(res["classes"], {"exact": 3, "over": 2, "under": 1, "no_is": 1})
        per = {r["is_name"]: r for r in csv.DictReader(open(os.path.join(out, "per_is.tsv")), delimiter="\t")}
        self.assertEqual(per["ISX"]["modal_extra_bp"], "9")
        # a second run into the same directory, and a wrong pin, are refused
        self.assertNotEqual(subprocess.run(cmd, capture_output=True).returncode, 0)
        cmd[cmd.index("--out") + 1] = out + "2"; cmd[cmd.index("--isfinder-sha256") + 1] = "0" * 64
        self.assertNotEqual(subprocess.run(cmd, capture_output=True).returncode, 0)
        self.assertFalse(os.path.exists(out + "2"))


if __name__ == "__main__":
    unittest.main()
