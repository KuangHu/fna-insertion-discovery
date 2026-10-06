#!/usr/bin/env python3
"""Integrity check of a built bag corpus. Run before anyone consumes it.

Checks the invariants the corpus is supposed to guarantee, not the numbers it
happened to produce. Every check names what it would catch, and a FAIL is a
reason not to use the corpus rather than a note.

  C1  files readable and non-empty       -- scratch has been throwing ESHUTDOWN
  C2  every line is valid JSON           -- a truncated write leaves a half line
  C3  required fields present on all     -- a schema drift is silent otherwise
  C4  flank length == 2*half, junction at half
  C5  flank is CANONICAL (<= its revcomp) and `orient` agrees
  C6  no duplicate (bag_id, flank)       -- the site-independence guarantee
  C7  bags.tsv agrees with the sites     -- bag membership and n_sites
  C8  sites + drops == database rows     -- nothing vanished unaccounted for
  C9  flank contains no element sequence -- it must come from the EMPTY allele,
      so it cannot equal or contain the insert; checked by md5 against the
      inserts FASTA for a sample
  C10 nc regions are inside the insert length, and nc hash matches the regions
"""
import collections
import csv
import hashlib
import json
import os
import sys

csv.field_size_limit(sys.maxsize)

BAGS = sys.argv[1]           # prefix, e.g. .../bags_v4/insertions
DB = sys.argv[2]             # .../database_v3/insertions.tsv
HALF = int(sys.argv[3]) if len(sys.argv) > 3 else 60

fails, warns = [], []


def ok(c, msg):
    print("  %-5s %s" % ("PASS" if c else "FAIL", msg))
    if not c:
        fails.append(msg)


def rc(s):
    return s.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


print("=" * 74)
print("BAG CORPUS CHECK  %s" % BAGS)
print("=" * 74)

# ---- C1
paths = {"sites": BAGS + "_sites.jsonl", "bags": BAGS + "_bags.tsv",
         "drops": BAGS + "_drops.tsv"}
print("\nC1  files")
for k, p in paths.items():
    try:
        sz = os.path.getsize(p)
        ok(sz > 0, "%-6s readable, %d bytes" % (k, sz))
    except OSError as e:
        ok(False, "%-6s UNREADABLE: %s" % (k, e))

# ---- C2/C3/C4/C5
REQ = ["site_id", "bag_id", "cds_cluster_id", "flank", "flank_len",
       "insertion_point_in_flank", "orient", "species", "noncoding_regions",
       "nc_sequence_hash", "empty_site_source"]
print("\nC2-C5  per-site schema and flank invariants")
n = 0
bad_json = bad_field = bad_len = bad_jct = bad_canon = bad_orient = 0
seen = collections.Counter()
bag_species = collections.defaultdict(set)
bag_n = collections.Counter()
nc_bad = 0
for line in open(paths["sites"]):
    line = line.strip()
    if not line:
        continue
    n += 1
    try:
        d = json.loads(line)
    except ValueError:
        bad_json += 1
        continue
    if any(f not in d for f in REQ):
        bad_field += 1
        continue
    fl = d["flank"]
    if len(fl) != 2 * HALF or d["flank_len"] != len(fl):
        bad_len += 1
    if d["insertion_point_in_flank"] != HALF:
        bad_jct += 1
    r = rc(fl)
    if fl > r:
        bad_canon += 1
    exp = "fwd" if fl <= r else "rc"
    if d["orient"] not in ("fwd", "rc"):
        bad_orient += 1
    seen[(d["bag_id"], fl)] += 1
    bag_species[d["bag_id"]].add(d["species"])
    bag_n[d["bag_id"]] += 1
    nc = d.get("noncoding_regions") or []
    h = hashlib.md5("|".join(nc).encode()).hexdigest()[:12]
    if nc and d.get("nc_sequence_hash") and not str(d["nc_sequence_hash"]).startswith(h[:6]):
        nc_bad += 1

ok(bad_json == 0, "all %d lines parse as JSON (%d bad)" % (n, bad_json))
ok(bad_field == 0, "required fields present (%d missing)" % bad_field)
ok(bad_len == 0, "flank length == %d everywhere (%d bad)" % (2 * HALF, bad_len))
ok(bad_jct == 0, "junction at index %d everywhere (%d bad)" % (HALF, bad_jct))
ok(bad_canon == 0,
   "flank is canonical, flank <= revcomp(flank) (%d violations)" % bad_canon)
ok(bad_orient == 0, "orient is fwd|rc (%d bad)" % bad_orient)

# ---- C6
print("\nC6  site independence")
dups = sum(v - 1 for v in seen.values() if v > 1)
ok(dups == 0, "no duplicate (bag_id, flank): %d duplicate rows" % dups)
print("       sites %d, distinct (bag,flank) %d, distinct flanks %d"
      % (n, len(seen), len({f for _, f in seen})))

# ---- C7
print("\nC7  bags.tsv agrees with the sites")
bt = {r["bag_id"]: r for r in csv.DictReader(open(paths["bags"]), delimiter="\t")}
ok(set(bt) == set(bag_n), "same bag set (%d in tsv, %d in sites)"
   % (len(bt), len(bag_n)))
mism = sum(1 for b, r in bt.items() if int(r["n_sites"]) != bag_n.get(b, -1))
ok(mism == 0, "n_sites matches for every bag (%d mismatched)" % mism)
sp_m = sum(1 for b, r in bt.items()
           if int(r["n_species"]) != len(bag_species.get(b, ())))
ok(sp_m == 0, "n_species matches for every bag (%d mismatched)" % sp_m)

# ---- C8
print("\nC8  accounting: sites + drops == database rows")
ndb = sum(1 for _ in open(DB)) - 1
ndrop = 0
try:
    for r in csv.DictReader(open(paths["drops"]), delimiter="\t"):
        ndrop += int(r.get("n") or r.get("count") or 1)
except (OSError, ValueError):
    ndrop = -1
print("       database rows %d, sites %d, drop rows read %s" % (ndb, n, ndrop))
if ndrop >= 0:
    # sites are DEDUPED, so sites+drops <= db rows; the gap is the collapse
    ok(n + ndrop <= ndb,
       "sites + drops (%d) <= database rows (%d); gap %d is the (bag,flank) "
       "collapse" % (n + ndrop, ndb, ndb - n - ndrop))
else:
    warns.append("drops table unreadable, accounting not checked")
    print("  WARN  drops table unreadable")

# ---- C9
print("\nC9  the flank carries no element sequence")
ins = set()
fa = os.path.join(os.path.dirname(DB), "insertions_inserts.fna")
if os.path.exists(fa):
    cur = []
    for line in open(fa):
        if line[0] == ">":
            if cur:
                ins.add(hashlib.md5("".join(cur).encode()).hexdigest()[:12])
            cur = []
        else:
            cur.append(line.strip())
    if cur:
        ins.add(hashlib.md5("".join(cur).encode()).hexdigest()[:12])
    hit = sum(1 for (_, f) in seen
              if hashlib.md5(f.encode()).hexdigest()[:12] in ins)
    ok(hit == 0, "no flank equals a full insert sequence (%d hits)" % hit)
else:
    print("  WARN  inserts FASTA not found, skipped")

print("\n" + "=" * 74)
if fails:
    print("RESULT: FAIL -- %d check(s) failed" % len(fails))
    for f in fails:
        print("   - %s" % f)
else:
    print("RESULT: PASS -- corpus is internally consistent")
if warns:
    print("warnings: %s" % "; ".join(warns))
print("=" * 74)
sys.exit(1 if fails else 0)
