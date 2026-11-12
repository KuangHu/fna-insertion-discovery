#!/usr/bin/env python3
"""Take six-frame rescue candidates back to the carrier contig and ask what they are.

A domain found inside a bare extraction is a fact about the extraction.  The
question it raises -- small CDS, or the amputated end of a large one? -- cannot
be answered inside the extraction, because the extraction is exactly what cut
it.  So every candidate is re-examined IN SITU, on the real carrier contig,
with the true upstream and downstream sequence continuing past the insertion
boundary.

WHAT IS DONE, IN ORDER
 1. The domain's interval is carried from the deduplicated representative back
    to one concrete occurrence.  The representative may be the reverse
    complement of that occurrence (see 220), so the transform is applied
    explicitly and the reading strand flips with it:
        orient fwd:  member[i] == rep[i]
        orient rev:  member[i] == complement(rep[L-1-i]);  '+' <-> '-'
    Getting this wrong would put the ORF search in the wrong frame on the
    wrong strand and quietly invent genes, so the mapped interval is VERIFIED
    by re-extracting it from the contig and checking it matches the
    representative slice base for base.  A mismatch is a hard failure for that
    candidate, never a warning.
 2. In that frame, walk outward on the contig to the nearest in-frame stop on
    each side.  The walk uses the contig, so it crosses the insertion boundary
    whenever the real ORF does.
 3. Inside the stop-to-stop span, collect in-frame start codons (table 11:
    ATG/GTG/TTG).  The start reported is the FIRST one after the upstream stop
    -- and it is reported as a CANDIDATE, with the full list of alternatives
    and the distance to the domain, never as a decision.  The brief is explicit
    that a start may not be set from the longest ORF alone, and a Pfam domain
    edge is not a CDS edge.
 4. Prodigal is re-run in situ, on the ORF plus 1 kb of real flank on each
    side, under the same meta/table-11 settings.  A gene it calls there that it
    did not call in the bare extraction is the strongest single piece of
    evidence that the extraction boundary, not the sequence, was the problem.
 5. The relationship to the insertion is RECORDED AND LEFT ALONE.  If the ORF
    runs past the insertion boundary, that fact is written down
    (`overlap_class`, `outside_5p_nt`, `outside_3p_nt`).  The insertion
    boundary is NOT moved to accommodate a gene, and an ORF whose coding
    majority lies outside the insertion is marked `likely_host_gene` rather
    than being credited to the element.

STATUS, assigned from evidence and nothing else:
   complete_candidate    stops on both sides, a start between them, nothing
                         running off a contig end
   start_uncertain       bounded, but no in-frame start after the upstream stop
   partial               the span reaches a contig end, so an end is unknown
   frameshift_or_conflict  in-situ Prodigal covers the domain in another frame
"""
import argparse, collections, json, os, sys

import pyrodigal
from Bio.Data import CodonTable

COMP = bytes.maketrans(b"ACGTacgt", b"TGCAtgca")
WD = "/global/scratch/users/kh36969/fna_ins_discovery"
STOPS = {"TAA", "TAG", "TGA"}
STARTS = {"ATG", "GTG", "TTG"}
FWD = dict(CodonTable.unambiguous_dna_by_id[11].forward_table)


def revcomp(s):
    return s.translate(COMP)[::-1]


def read_tsv(p):
    with open(p) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        for ln in f:
            yield dict(zip(hdr, ln.rstrip("\n").split("\t")))


def load_contigs(path, wanted):
    out, rid, parts = {}, None, []
    with open(path, "rb") as f:
        for line in f:
            if line[:1] == b">":
                if rid in wanted:
                    out[rid] = b"".join(parts).upper()
                rid = line[1:].split()[0].decode()
                parts = []
            else:
                parts.append(line.strip())
    if rid in wanted:
        out[rid] = b"".join(parts).upper()
    return out


def walk_stop(ctg, pos, step, cap):
    """From in-frame codon start `pos`, step by +/-3 to the first stop codon.

    Returns (codon_start, 'stop') or (edge_position, 'contig_end'/'cap')."""
    n, p = 0, pos
    while True:
        p += step
        n += 1
        if p < 0 or p + 3 > len(ctg):
            return p - step, "contig_end"
        if n * 3 > cap:
            return p, "cap"
        if ctg[p:p + 3].decode() in STOPS:
            return p, "stop"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cand", required=True)
    ap.add_argument("--members", required=True)
    ap.add_argument("--browse", required=True)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--repfa", required=True, help="pilot_sample.fna (representatives)")
    ap.add_argument("--wd", default=WD)
    ap.add_argument("--out", required=True, help="output prefix")
    ap.add_argument("--cap", type=int, default=30000)
    ap.add_argument("--flank", type=int, default=1000)
    a = ap.parse_args()

    # representative sequences, for the mapping check
    rep = {}
    rid, parts = None, []
    with open(a.repfa) as f:
        for line in f:
            if line[0] == ">":
                if rid:
                    rep[rid] = "".join(parts).upper()
                t = line[1:].split()[0]
                rid = t.split("|", 1)[1] if "|" in t else t
                parts = []
            else:
                parts.append(line.strip())
    if rid:
        rep[rid] = "".join(parts).upper()

    cands = list(read_tsv(a.cand))
    keys = {c["seq_key"] for c in cands}
    print(f"rescue candidates {len(cands)} over {len(keys)} sequences", flush=True)

    # one occurrence per key: smallest occurrence_id, matching how 220 chose the rep
    pick = {}
    for r in read_tsv(a.members):
        if r["seq_key"] not in keys:
            continue
        cur = pick.get(r["seq_key"])
        if cur is None or r["occurrence_id"] < cur["occurrence_id"]:
            pick[r["seq_key"]] = r

    need_occ = {v["occurrence_id"]: k for k, v in pick.items()}
    occ = {}
    for r in read_tsv(a.browse):
        if r["occurrence_id"] in need_occ:
            occ[r["occurrence_id"]] = r
    mani = {}
    with open(a.manifest) as f:
        for ln in f:
            j = json.loads(ln)
            mani[(j["assembly_accession_version"], j["fasta_record_id"])] = j["fna_path"]

    # group by genome file so each assembly is read once
    by_file = collections.defaultdict(list)
    for c in cands:
        k = c["seq_key"]
        p = pick.get(k)
        if p is None:
            print(f"  no member for {k}", flush=True)
            continue
        o = occ[p["occurrence_id"]]
        fp = mani.get((o["carrier_assembly_accession_version"], o["carrier_fasta_record_id"]))
        if fp is None:
            print(f"  no manifest entry for {p['occurrence_id']}", flush=True)
            continue
        by_file[fp].append((c, p, o))

    gf = pyrodigal.GeneFinder(
        meta=True,
        metagenomic_bins=pyrodigal.MetagenomicBins(
            [b for b in pyrodigal.METAGENOMIC_BINS if b.training_info.translation_table == 11]),
        closed=False, mask=False)

    COLS = ["seq_key", "occurrence_id", "pfam_acc", "pfam_name", "assembly", "contig",
            "contig_len", "insert_start0", "insert_end0", "orient_to_rep",
            "domain_contig_start0", "domain_contig_end0", "read_strand",
            "up_stop_kind", "dn_stop_kind", "orf_start0", "orf_end0", "orf_len_nt",
            "n_starts_in_orf", "first_start0", "start_codon", "dist_start_to_domain_nt",
            "cds_start0", "cds_end0", "cds_len_aa", "overlap_class",
            "outside_5p_nt", "outside_3p_nt", "frac_inside_insertion",
            "insitu_prodigal", "insitu_same_frame", "status", "likely_host_gene",
            "map_check"]
    out = open(a.out + "_rescue.tsv", "w")
    out.write("\t".join(COLS) + "\n")
    faa = open(a.out + "_rescued.faa", "w")
    cnt = collections.Counter()

    for fp, items in sorted(by_file.items()):
        full = os.path.join(a.wd, fp)
        wanted = {o["carrier_fasta_record_id"] for (_, _, o) in items}
        try:
            ctgs = load_contigs(full, wanted)
        except OSError as e:
            print(f"  unreadable {fp}: {e}", flush=True)
            cnt["genome_unreadable"] += len(items)
            continue
        for (c, p, o) in items:
            ctg = ctgs.get(o["carrier_fasta_record_id"])
            if ctg is None:
                cnt["contig_missing"] += 1
                continue
            CL = len(ctg)
            ins0, ins1 = int(o["insert_start0"]), int(o["insert_end0"])
            key = c["seq_key"]
            R = rep[key]
            Lrep = len(R)
            d0, d1 = int(c["dna_start0"]), int(c["dna_end0"])
            strand = c["strand"]

            # representative -> member -> contig
            if p["orient_to_rep"] == "fwd":
                m0, m1, rstrand = d0, d1, strand
            else:
                m0, m1 = Lrep - d1, Lrep - d0
                rstrand = "-" if strand == "+" else "+"
            c0, c1 = ins0 + m0, ins0 + m1

            # the mapping is verified against the DNA, not assumed
            got = ctg[c0:c1].decode()
            exp = R[d0:d1]
            if p["orient_to_rep"] == "rev":
                exp = revcomp(exp.encode()).decode()
            mapchk = "ok" if got == exp else "MISMATCH"
            if mapchk != "ok":
                cnt["map_mismatch"] += 1
                out.write("\t".join([key, p["occurrence_id"], c["pfam_acc"], c["pfam_name"],
                                     o["carrier_assembly_accession_version"],
                                     o["carrier_fasta_record_id"], str(CL), str(ins0), str(ins1),
                                     p["orient_to_rep"], str(c0), str(c1), rstrand] +
                                    [""] * (len(COLS) - 14) + ["MISMATCH"]) + "\n")
                continue

            # ORF walk, in the domain's own frame, on the contig
            work = ctg if rstrand == "+" else revcomp(ctg)
            if rstrand == "+":
                w0, w1 = c0, c1
            else:
                w0, w1 = CL - c1, CL - c0
            up_pos, up_kind = walk_stop(work, w0, -3, a.cap)
            dn_pos, dn_kind = walk_stop(work, w1 - 3, +3, a.cap)
            orf0 = up_pos + 3 if up_kind == "stop" else up_pos
            orf1 = dn_pos if dn_kind == "stop" else min(len(work), dn_pos + 3)
            starts = [i for i in range(orf0, w0 + 1, 3)
                      if work[i:i + 3].decode() in STARTS]
            first = starts[0] if starts else -1

            if first >= 0:
                cds0, cds1 = first, (dn_pos + 3 if dn_kind == "stop" else orf1)
                prot = []
                for i in range(cds0, cds1 - 2, 3):
                    cod = work[i:i + 3].decode()
                    if cod in STOPS:
                        break
                    prot.append(FWD.get(cod, "X"))
                if prot:
                    prot[0] = "M"
                prot = "".join(prot)
            else:
                cds0 = cds1 = -1
                prot = ""

            # back to contig forward coordinates for reporting
            if rstrand == "+":
                f_orf0, f_orf1 = orf0, orf1
                f_cds0, f_cds1 = cds0, cds1
                f_first = first
            else:
                f_orf0, f_orf1 = len(work) - orf1, len(work) - orf0
                f_cds0, f_cds1 = (len(work) - cds1, len(work) - cds0) if first >= 0 else (-1, -1)
                f_first = len(work) - first - 3 if first >= 0 else -1

            lo, hi = (f_cds0, f_cds1) if first >= 0 else (f_orf0, f_orf1)
            inside = max(0, min(hi, ins1) - max(lo, ins0))
            span = max(1, hi - lo)
            out5 = max(0, ins0 - lo)
            out3 = max(0, hi - ins1)
            if out5 and out3:
                ocls = "spans_both_boundaries"
            elif out5:
                ocls = "spans_5prime_boundary"
            elif out3:
                ocls = "spans_3prime_boundary"
            else:
                ocls = "inside_insertion"
            host = int(inside / span < 0.5)

            # in-situ Prodigal, with real flanks
            s = max(0, lo - a.flank)
            e = min(CL, hi + a.flank)
            win = ctg[s:e].decode()
            insitu, same_frame = "none", 0
            if len(win) >= 60:
                for g in gf.find_genes(win):
                    gs, ge = s + g.begin - 1, s + g.end
                    if ge <= c0 or gs >= c1:
                        continue
                    gstr = "+" if g.strand > 0 else "-"
                    if gstr == rstrand:
                        gf0 = gs if gstr == "+" else CL - ge
                        if gf0 % 3 == (w0 % 3):
                            insitu, same_frame = "covers_same_frame", 1
                            break
                    insitu = "covers_other_frame"

            if insitu == "covers_other_frame" and not same_frame:
                status = "frameshift_or_conflict"
            elif up_kind != "stop" or dn_kind != "stop":
                status = "partial"
            elif first < 0:
                status = "start_uncertain"
            else:
                status = "complete_candidate"
            cnt[status] += 1
            cnt["insitu_" + insitu] += 1

            pid = f"R|{key}|{c['pfam_acc']}|{c0}-{c1}"
            if prot:
                faa.write(f">{pid} {o['carrier_assembly_accession_version']}"
                          f"|{o['carrier_fasta_record_id']}|{f_cds0}-{f_cds1}"
                          f"|{rstrand} status={status}\n{prot}\n")
            out.write("\t".join(str(x) for x in [
                key, p["occurrence_id"], c["pfam_acc"], c["pfam_name"],
                o["carrier_assembly_accession_version"], o["carrier_fasta_record_id"],
                CL, ins0, ins1, p["orient_to_rep"], c0, c1, rstrand,
                up_kind, dn_kind, f_orf0, f_orf1, f_orf1 - f_orf0,
                len(starts), f_first, work[first:first + 3].decode() if first >= 0 else "",
                (w0 - first) if first >= 0 else "",
                f_cds0, f_cds1, len(prot), ocls, out5, out3,
                f"{inside/span:.3f}", insitu, same_frame, status, host, mapchk]) + "\n")

    out.close()
    faa.close()
    with open(a.out + "_rescue.json", "w") as f:
        json.dump(dict(counts=dict(cnt), candidates=len(cands),
                       sequences=len(keys), cap_nt=a.cap, flank_nt=a.flank), f, indent=2)
    print(json.dumps(dict(cnt), indent=2))


if __name__ == "__main__":
    main()
