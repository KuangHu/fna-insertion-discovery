#!/usr/bin/env python3
"""Put the two routes' Pfam hits on the same DNA ruler and classify the difference.

A Pfam hit is scored on a PROTEIN.  Comparing route A and route B means first
putting both back onto the DNA they came from, because that is the only
coordinate system the two share.

   strand '+' :  dna_start0 = nt_start0 + (ali_from-1)*3
                 dna_end0   = nt_start0 +  ali_to    *3
   strand '-' :  dna_end0   = nt_end0   - (ali_from-1)*3
                 dna_start0 = nt_end0   -  ali_to    *3

FRAME IS PART OF THE IDENTITY OF A HIT, NOT A DETAIL.  Two intervals can
overlap on the DNA and still be different coding claims.  The frame key used
here is (strand, offset), offset measured along the strand's OWN sequence:
   '+' -> nt_start0 % 3          '-' -> (L - nt_end0) % 3
Classifying a '-' hit against a '+' CDS as "the same gene" because the
intervals overlap would be the same mistake as mixing one frame's left end
with another frame's right end.

THE FIVE OUTCOMES ARE THE BRIEF'S, UNCHANGED:
  within_existing_cds           inside a Prodigal CDS, same strand and frame
                                -> annotation merged, NOT a new CDS
  cds_boundary_extension        same frame, but the domain runs past the
                                predicted boundary -> boundary is the suspect
  sixframe_rescue_candidate     no Prodigal CDS at that place at all
  frame_conflict                overlaps a CDS in a different frame or strand
                                -> conflict, never silently a new CDS
  (orthogonal) support_level    full_domain  hmm_cov >= 0.80
                                domain_fragment  otherwise
                                A fragment is kept as EVIDENCE OF CODING, not
                                promoted to a CDS.  This is a second axis, not
                                a sixth category, because a fragment can occur
                                in any of the four geometries.

WHAT THIS SCRIPT DOES NOT DO.  It does not call a start codon, does not move
any boundary, and does not create a CDS.  A rescue candidate leaves here as a
question for 225, which goes back to the carrier contig and the real flanks.
A Pfam domain edge is NOT a CDS edge.
"""
import argparse, collections, json, os, sys

COV_FULL = 0.80


def read_tsv(p):
    with open(p) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        for ln in f:
            yield dict(zip(hdr, ln.rstrip("\n").split("\t")))


def frame_key(strand, s0, e0, L):
    return (strand, s0 % 3 if strand == "+" else (L - e0) % 3)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--calls", required=True)
    ap.add_argument("--pfam", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--out", required=True, help="output prefix")
    ap.add_argument("--cov-full", type=float, default=COV_FULL)
    a = ap.parse_args()

    # sequence lengths are needed for the '-' strand frame key
    L = {}
    rid = None
    n = 0
    with open(a.fasta) as f:
        for line in f:
            if line[0] == ">":
                rid = line[1:].split()[0]
                rid = rid.split("|", 1)[1] if "|" in rid else rid
                L[rid] = 0
            else:
                L[rid] += len(line.strip())

    calls = {}
    cds_by_seq = collections.defaultdict(list)
    for r in read_tsv(a.calls):
        calls[r["call_id"]] = r
        if r["route"] == "prodigal":
            s0, e0 = int(r["nt_start0"]), int(r["nt_end0"])
            cds_by_seq[r["seq_key"]].append(
                (s0, e0, r["strand"], frame_key(r["strand"], s0, e0, L[r["seq_key"]]),
                 r["call_id"], int(r["partial_begin"]), int(r["partial_end"])))

    rows, cnt = [], collections.Counter()
    fam_by_seq = collections.defaultdict(lambda: collections.defaultdict(set))
    for h in read_tsv(a.pfam):
        c = calls.get(h["protein_id"])
        if c is None:
            sys.exit(f"FATAL: pfam hit on unknown call {h['protein_id']}")
        key, strand = c["seq_key"], c["strand"]
        s0, e0 = int(c["nt_start0"]), int(c["nt_end0"])
        af, at = int(h["ali_from"]), int(h["ali_to"])
        if strand == "+":
            d0, d1 = s0 + (af - 1) * 3, s0 + at * 3
        else:
            d1, d0 = e0 - (af - 1) * 3, e0 - at * 3
        fk = frame_key(strand, s0, e0, L[key])
        cov = float(h["hmm_cov"])
        support = "full_domain" if cov >= a.cov_full else "domain_fragment"
        fam_by_seq[key][h["route"]].add(h["pfam_acc"].split(".")[0])

        geom, partner, detail = "", "", ""
        if h["route"] == "prodigal":
            geom, partner = "prodigal_cds", c["call_id"]
        else:
            same, other = [], []
            for (cs, ce, cst, cfk, cid, pb, pe) in cds_by_seq[key]:
                if ce <= d0 or cs >= d1:
                    continue
                (same if cfk == fk else other).append((cs, ce, cid, pb, pe))
            if same:
                cs, ce, cid, pb, pe = same[0]
                partner = cid
                if d0 >= cs and d1 <= ce:
                    geom = "within_existing_cds"
                else:
                    geom = "cds_boundary_extension"
                    detail = (f"domain {d0}-{d1} vs cds {cs}-{ce}; "
                              f"left_out={max(0, cs-d0)} right_out={max(0, d1-ce)}; "
                              f"cds_partial={pb}{pe}")
            elif other:
                geom = "frame_conflict"
                partner = ",".join(x[2] for x in other)
            else:
                geom = "sixframe_rescue_candidate"
        cnt[(h["route"], geom, support)] += 1
        rows.append([h["protein_id"], h["route"], key, h["pfam_acc"], h["pfam_name"],
                     strand, fk[1], d0, d1, h["bitscore"], h["i_evalue"],
                     h["hmm_cov"], support, geom, partner, detail])

    with open(a.out + "_hits_on_dna.tsv", "w") as f:
        f.write("protein_id\troute\tseq_key\tpfam_acc\tpfam_name\tstrand\tframe\t"
                "dna_start0\tdna_end0\tbitscore\ti_evalue\thmm_cov\tsupport_level\t"
                "geometry\tpartner_cds\tdetail\n")
        for r in sorted(rows, key=lambda r: (r[2], r[7], r[3])):
            f.write("\t".join(str(x) for x in r) + "\n")

    with open(a.out + "_rescue_candidates.tsv", "w") as f:
        f.write("seq_key\tpfam_acc\tpfam_name\tstrand\tframe\tdna_start0\tdna_end0\t"
                "hmm_cov\tsupport_level\tbitscore\tprotein_id\n")
        for r in sorted(rows, key=lambda r: (r[2], r[7])):
            if r[13] == "sixframe_rescue_candidate":
                f.write("\t".join(str(x) for x in
                                  [r[2], r[3], r[4], r[5], r[6], r[7], r[8],
                                   r[11], r[12], r[9], r[0]]) + "\n")

    # family-level agreement per insertion sequence
    agree = collections.Counter()
    for key, d in fam_by_seq.items():
        A, B = d.get("prodigal", set()), d.get("sixframe", set())
        agree["both_families"] += len(A & B)
        agree["prodigal_only_families"] += len(A - B)
        agree["sixframe_only_families"] += len(B - A)
        if A and B:
            agree["seqs_with_both_routes"] += 1
        elif A:
            agree["seqs_prodigal_only"] += 1
        elif B:
            agree["seqs_sixframe_only"] += 1

    summ = {"geometry_x_support": {f"{r}|{g}|{s}": v for (r, g, s), v in sorted(cnt.items())},
            "family_agreement": dict(agree),
            "sequences_with_any_hit": len(fam_by_seq),
            "cov_full_threshold": a.cov_full}
    with open(a.out + "_compare.json", "w") as f:
        json.dump(summ, f, indent=2)
    print(json.dumps(summ, indent=2))


if __name__ == "__main__":
    main()
