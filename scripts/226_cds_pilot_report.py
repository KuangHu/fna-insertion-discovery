#!/usr/bin/env python3
"""Round-1 report: what the two routes agree on, what six-frame adds, what it costs.

Every headline number here is paired with the control that makes it readable.
"Six-frame added N domains" means nothing until the same pipeline has been run
on composition-matched shuffles of the same sequences; "Prodigal missed M
genes" means nothing until real genes from the same organisms have been cut
the same way and fed back in.

POSITIVE-CONTROL RECOVERY IS JUDGED AGAINST THE STORED TRUE PROTEIN, never
against the route's own self-consistency.  The expected protein of a truncated
control is re-derived from the control's DNA in the gene's own frame
    frame offset in the extraction = (-trim5) mod 3
and checked to be a suffix of the true protein before it is used.  If that
check fails the control is reported as BROKEN rather than silently scored.

NOTHING IS DISCARDED FOR LACKING A PFAM HIT.  Prodigal CDS with no domain are
counted and kept; they are predictions, and Pfam's silence is Pfam's limit,
not evidence against the gene.
"""
import argparse, collections, json, os

from Bio.Data import CodonTable

FWD = dict(CodonTable.unambiguous_dna_by_id[11].forward_table)
STOPS = {"TAA", "TAG", "TGA"}


def read_tsv(p):
    if not os.path.exists(p):
        return
    with open(p) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        for ln in f:
            yield dict(zip(hdr, ln.rstrip("\n").split("\t")))


def fasta(p):
    rid, parts = None, []
    with open(p) as f:
        for line in f:
            if line[0] == ">":
                if rid:
                    yield rid, "".join(parts)
                rid = line[1:].split()[0]
                parts = []
            else:
                parts.append(line.strip())
    if rid:
        yield rid, "".join(parts)


def translate(dna, off):
    out = []
    for i in range(off, len(dna) - 2, 3):
        c = dna[i:i + 3]
        if c in STOPS:
            break
        out.append(FWD.get(c, "X"))
    return "".join(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--w", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--rescued-pfam", default="")
    a = ap.parse_args()
    W = a.w
    R = {}

    # ---------- pilot ---------------------------------------------------------
    calls = list(read_tsv(f"{W}/pilot_calls.tsv"))
    nseq = len({c["seq_key"] for c in calls})
    hits = list(read_tsv(f"{W}/pilot_hits_on_dna.tsv"))
    prot_with_hit = {h["protein_id"] for h in hits}
    prod_calls = [c for c in calls if c["route"] == "prodigal"]
    R["pilot"] = {
        "sequences": nseq,
        "prodigal_cds": len(prod_calls),
        "prodigal_cds_partial_either_end": sum(
            1 for c in prod_calls if c["partial_begin"] == "1" or c["partial_end"] == "1"),
        "prodigal_cds_with_pfam": len({c["call_id"] for c in prod_calls
                                       if c["call_id"] in prot_with_hit}),
        "prodigal_cds_without_pfam_RETAINED": len(
            [c for c in prod_calls if c["call_id"] not in prot_with_hit]),
        "sixframe_fragments_ge20aa": sum(1 for c in calls if c["route"] == "sixframe"),
        "sixframe_fragments_with_pfam": len({h["protein_id"] for h in hits
                                             if h["route"] == "sixframe"}),
    }
    geom = collections.Counter()
    for h in hits:
        if h["route"] == "sixframe":
            geom[(h["geometry"], h["support_level"])] += 1
    R["sixframe_hit_geometry"] = {f"{g}|{s}": n for (g, s), n in sorted(geom.items())}

    fam = collections.defaultdict(lambda: collections.defaultdict(set))
    for h in hits:
        fam[h["seq_key"]][h["route"]].add(h["pfam_acc"].split(".")[0])
    both = sum(len(d.get("prodigal", set()) & d.get("sixframe", set())) for d in fam.values())
    aonly = sum(len(d.get("prodigal", set()) - d.get("sixframe", set())) for d in fam.values())
    bonly = sum(len(d.get("sixframe", set()) - d.get("prodigal", set())) for d in fam.values())
    R["family_instances_per_sequence"] = {
        "supported_by_both_routes": both,
        "prodigal_only": aonly,
        "sixframe_only": bonly,
        "sequences_with_any_domain": len(fam),
    }

    # ---------- rescue --------------------------------------------------------
    resc = list(read_tsv(f"{W}/pilot_rescue.tsv"))
    if resc:
        st = collections.Counter(r["status"] for r in resc)
        R["context_rescue"] = {
            "candidates_examined": len(resc),
            "status": dict(st),
            "mapping_mismatch": sum(1 for r in resc if r["map_check"] != "ok"),
            "orf_crosses_insertion_boundary": sum(
                1 for r in resc if r.get("overlap_class", "") not in ("inside_insertion", "")),
            "likely_host_gene": sum(1 for r in resc if r.get("likely_host_gene") == "1"),
            "insitu_prodigal_recovers_same_frame": sum(
                1 for r in resc if r.get("insitu_same_frame") == "1"),
        }

    # ---------- did in-situ context actually COMPLETE the domain? -------------
    if a.rescued_pfam and os.path.exists(a.rescued_pfam):
        before = collections.defaultdict(float)
        for c in read_tsv(f"{W}/pilot_rescue_candidates.tsv"):
            k = (c["seq_key"], c["pfam_acc"].split(".")[0])
            before[k] = max(before[k], float(c["hmm_cov"]))
        after = collections.defaultdict(float)
        for h in read_tsv(a.rescued_pfam):
            pid = h["protein_id"].split("|")
            if len(pid) < 3:
                continue
            k = (pid[1], pid[2].split(".")[0])
            after[k] = max(after[k], float(h["hmm_cov"]))
        # A candidate with no rescued protein (no start found, so 225 emitted
        # nothing) is NOT a domain that was lost by rescuing -- it is a domain
        # that was never re-tested. Conflating the two would read as evidence
        # against context rescue when it is evidence of an unresolved start.
        produced = set()
        for h in read_tsv(f"{W}/pilot_rescue.tsv"):
            if int(h.get("cds_len_aa") or 0) > 0:
                produced.add((h["seq_key"], h["pfam_acc"].split(".")[0]))
        shared = [k for k in before if k in after]
        R["domain_completion_by_context"] = {
            "candidate_domains": len(before),
            "rescued_protein_rescanned_and_hit": len(shared),
            "coverage_increased": sum(1 for k in shared if after[k] > before[k] + 1e-9),
            "fragment_became_full_domain": sum(
                1 for k in shared if before[k] < 0.80 <= after[k]),
            "no_rescued_protein_to_test": sum(1 for k in before if k not in produced),
            "protein_produced_but_family_no_longer_hits": sum(
                1 for k in before if k in produced and k not in after),
            "median_cov_before": round(sorted(before[k] for k in shared)[len(shared)//2], 3)
            if shared else None,
            "median_cov_after": round(sorted(after[k] for k in shared)[len(shared)//2], 3)
            if shared else None,
        }

    # ---------- controls ------------------------------------------------------
    truth = {t["control_id"]: t for t in read_tsv(f"{W}/controls_truth.tsv")}
    cseq = dict(fasta(f"{W}/controls.fna")) if os.path.exists(f"{W}/controls.fna") else {}
    ccalls = collections.defaultdict(list)
    for c in read_tsv(f"{W}/controls_calls.tsv"):
        ccalls[c["seq_key"]].append(c)
    cprot = {}
    for p in (f"{W}/controls_prodigal.faa", f"{W}/controls_sixframe.faa"):
        if os.path.exists(p):
            cprot.update(dict(fasta(p)))
    chits = collections.defaultdict(set)
    for h in read_tsv(f"{W}/controls_pfam.tsv"):
        chits[h["seq_key"]].add((h["route"], h["pfam_acc"].split(".")[0]))

    neg = collections.Counter()
    pos = collections.defaultdict(lambda: collections.Counter())
    broken = 0
    for cid, t in truth.items():
        key = cid.split("|", 1)[1]
        routes_hit = {r for r, _ in chits.get(key, ())}
        if t["kind"].startswith("neg"):
            neg[t["kind"] + "_n"] += 1
            if "prodigal" in routes_hit:
                neg[t["kind"] + "_prodigal_hit"] += 1
            if "sixframe" in routes_hit:
                neg[t["kind"] + "_sixframe_hit"] += 1
            continue
        P = t["true_protein"]
        dna = cseq.get(cid, "")
        off = (-int(t["trim5_nt"])) % 3
        exp = translate(dna, off)
        if not P.endswith(exp[1:]) and exp[1:] not in P:
            broken += 1
            pos[t["kind"]]["BROKEN_expected_protein"] += 1
            continue
        k = t["kind"]
        pos[k]["n"] += 1
        for route in ("prodigal", "sixframe"):
            got = [cprot.get(c["call_id"], "") for c in ccalls.get(key, [])
                   if c["route"] == route]
            ok = any(exp in g or (len(g) >= 20 and g in exp) for g in got)
            pos[k][route + "_recovers_protein"] += int(ok)
            pos[k][route + "_any_call"] += int(bool(got))
        if "prodigal" in routes_hit:
            pos[k]["prodigal_pfam"] += 1
        if "sixframe" in routes_hit:
            pos[k]["sixframe_pfam"] += 1
    R["negative_controls"] = dict(neg)
    R["positive_controls"] = {k: dict(v) for k, v in sorted(pos.items())}
    R["positive_controls_broken_expectation"] = broken

    # false-discovery framing: per-sequence domain rate, real vs shuffled
    real_rate = len(fam) / max(1, nseq)
    nd = neg.get("neg_dinuc_sixframe_hit", 0) / max(1, neg.get("neg_dinuc_n", 1))
    nm = neg.get("neg_mono_sixframe_hit", 0) / max(1, neg.get("neg_mono_n", 1))
    R["per_sequence_domain_rate"] = {
        "real_insertions_any_route": round(real_rate, 4),
        "dinucleotide_shuffle_sixframe": round(nd, 4),
        "mononucleotide_shuffle_sixframe": round(nm, 4),
    }

    with open(a.out + ".json", "w") as f:
        json.dump(R, f, indent=2)
    print(json.dumps(R, indent=2))


if __name__ == "__main__":
    main()
