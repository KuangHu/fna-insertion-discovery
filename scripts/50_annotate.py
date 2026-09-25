#!/usr/bin/env python3
"""Stage 50 -- annotation, deliberately LAST, reported as SEPARATE EVIDENCE AXES.

Discovery in stages 20-40 never looks at a transposase model, so the catalogue
holds whatever is structurally mobile: known IS, unknown IS, prophage,
integrative elements, rrn operons, and the odd artefact. This stage only asks
WHAT each element is. It can never add or remove a candidate on family grounds.

The point of the rewrite: a single collapsed `class_call` threw away the
evidence that produced it, and an if/elif chain silently let one axis mask
another (an element with BOTH rRNA and a transposase domain reported only the
first branch tested). Every axis is now its own column, and the class is a
convenience label derived from them -- never a substitute for reading them.

EVIDENCE AXES
  structural  (carried in from the discovery stage, never recomputed here)
      n_copies, n_distinct_loci, boundary_dispersion, S_mobility, length
  ribosomal   barrnap                 -> is this an rrn repeat?
  coding      prodigal                -> longest ORF / length
  known IS    ISEScan's own ISfinder-derived pHMM library (clusters.faa.hmm)
              -> is this a KNOWN IS that the genome-level ISEScan run missed?
  domains     any HMM you pass with --hmm (Pfam-A, or the IS110 DEDD/Tnp20 pair)
              -> transposase / integrase / recombinase evidence

CLASSES (derived, in priority order, and every one keeps its evidence columns)
  ribosomal_repeat            rRNA genes present and ORF-poor
  known_IS_missed_by_ISEScan  hits the IS pHMM library though genome-level
                              ISEScan called nothing here -- NOT novel, but
                              direct proof that annotation-free discovery finds
                              known mobile elements a conventional caller missed
  known_other_MGE             integrase/recombinase/phage domain, not IS
  segmental_duplication       flank homology high / boundaries not sharp
  unknown_mobile_element_candidate
                              ORF-rich, no known hit, strong mobility evidence
  noncanonical_mobile_candidate
                              ORF-POOR but mobility evidence strong. Deliberately
                              NOT called an artefact: the whole point of the
                              project is elements that are not conventional
                              transposases. This class is a promotion, not a bin.
  unclassified                everything else

TSD/TIR are never inputs here. IS110/IS1111 makes neither.
"""
import argparse
import os
import re
import subprocess
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib.util import log, run                                   # noqa: E402


def read_fasta_lengths(path):
    lens, name = {}, None
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].split()[0]
            lens[name] = 0
        elif name:
            lens[name] += len(line.strip())
    return lens


def barrnap_hits(fa, outdir, barrnap, threads):
    gff = os.path.join(outdir, "barrnap.gff")
    try:
        with open(gff, "w") as fh:
            run([barrnap, "--threads", str(threads), "--quiet", fa],
                stdout=fh, capture_output=False)
    except Exception as exc:
        log("barrnap unavailable (%s) -- rRNA axis is UNKNOWN, not negative" % exc)
        return None
    hits = defaultdict(list)
    for line in open(gff):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) >= 9:
            m = re.search(r"Name=([^;]+)", f[8])
            hits[f[0]].append(m.group(1) if m else "rRNA")
    return hits


def prodigal_proteins(fa, outdir, prodigal):
    faa = os.path.join(outdir, "elements.faa")
    run([prodigal, "-i", fa, "-a", faa, "-p", "meta", "-q",
         "-o", os.path.join(outdir, "prodigal.gbk")])
    # Take the ORF length from prodigal's own coordinates:
    #   >elem_3 # 1610 # 3730 # 1 # ID=...
    # Summing sequence lines instead is a trap -- prodigal wraps protein FASTA
    # at 60 aa, so a per-line max silently reports 180 bp for every element.
    lens = defaultdict(int)
    for line in open(faa):
        if not line.startswith(">"):
            continue
        f = line[1:].split("#")
        if len(f) < 3:
            continue
        name = f[0].strip().rsplit("_", 1)[0]
        try:
            lens[name] = max(lens[name], abs(int(f[2]) - int(f[1])) + 1)
        except ValueError:
            continue
    return faa, lens


def hmm_hits(faa, hmm, tbl, hmmsearch, threads, evalue, cut_ga=False):
    """protein -> {(profile, evalue)}. Element id is the prodigal name minus _N."""
    cmd = [hmmsearch, "--cpu", str(threads), "--tblout", tbl]
    cmd += ["--cut_ga"] if cut_ga else ["-E", str(evalue)]
    cmd += [hmm, faa]
    try:
        run(cmd, capture_output=True)
    except Exception as exc:
        log("hmmsearch failed on %s (%s)" % (os.path.basename(hmm), exc))
        return None
    hits = defaultdict(list)
    for line in open(tbl):
        if line.startswith("#"):
            continue
        f = line.split()
        if len(f) > 4:
            hits[f[0].rsplit("_", 1)[0]].append((f[2], float(f[4])))
    return hits


def best(hitlist, n=3):
    if not hitlist:
        return ".", ""
    s = sorted(set(hitlist), key=lambda t: t[1])
    return ",".join(p for p, _ in s[:n]), "%.1e" % s[0][1]


MGE_WORDS = re.compile(
    r"integrase|recombinase|resolvase|transposase|phage|capsid|terminase|"
    r"portal|tail|relaxase|mob[A-Z]?|conjug|excision|invertase|tyrosine_recomb|"
    r"serine_recomb|rve|zinc_ribbon|rusa|ning|holliday|antirepressor|"
    r"gp5_trimer|tape_meas|nu1", re.I)
IS_WORDS = re.compile(r"transposase|tnp|^is\d|dde|dedd", re.I)
# A confident Pfam identity that is NOT an MGE domain still IDENTIFIES the
# element -- it is not a novel candidate. The dominant case in the first real
# run was Rhs / YD-repeat polymorphic toxin systems (RHS_repeat, TEN_YD-shell,
# DUF6531): genuinely repeated, genuinely not a transposon. Letting those fall
# through to "unknown_mobile_element_candidate" would have inflated the
# discovery pool with a well-characterised family. A DUF is not an identity, so
# it never rescues an element on its own.
def informative(dom):
    return dom != "." and any(not d.upper().startswith("DUF")
                              for d in dom.split(","))


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--elements", required=True, help="element FASTA")
    ap.add_argument("--structural", default=None,
                    help="TSV of structural evidence keyed by the FASTA id "
                         "(e.g. bench60_novel_pool.tsv or universal_elements.tsv)")
    ap.add_argument("--structural-key", default="family_id")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--is-hmm", default=None,
                    help="ISEScan's clusters.faa.hmm -- the known-IS axis")
    ap.add_argument("--hmm", default=None,
                    help="extra HMM db (Pfam-A, or the IS110 DEDD/Tnp20 pair)")
    ap.add_argument("--hmm-cut-ga", action="store_true",
                    help="use Pfam gathering thresholds instead of an E-value")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--evalue", type=float, default=1e-5)
    ap.add_argument("--min-orf-frac", type=float, default=0.4)
    ap.add_argument("--strong-s", type=float, default=0.70,
                    help="S_mobility at/above this counts as strong mobility")
    ap.add_argument("--barrnap", default="barrnap")
    ap.add_argument("--prodigal", default="prodigal")
    ap.add_argument("--hmmsearch", default="hmmsearch")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    lens = read_fasta_lengths(args.elements)
    struct = {}
    if args.structural and os.path.exists(args.structural):
        with open(args.structural) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            ki = hdr.index(args.structural_key)
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) == len(hdr):
                    struct[f[ki]] = dict(zip(hdr, f))

    log("%d elements" % len(lens))
    rrna = barrnap_hits(args.elements, args.outdir, args.barrnap, args.threads)
    faa, orf = prodigal_proteins(args.elements, args.outdir, args.prodigal)
    ishits = hmm_hits(faa, args.is_hmm, os.path.join(args.outdir, "is_hmm.tbl"),
                      args.hmmsearch, args.threads, args.evalue) \
        if args.is_hmm else None
    dhits = hmm_hits(faa, args.hmm, os.path.join(args.outdir, "dom_hmm.tbl"),
                     args.hmmsearch, args.threads, args.evalue,
                     cut_ga=args.hmm_cut_ga) if args.hmm else None

    cols = ["element_id", "length", "n_copies", "n_distinct_loci",
            "boundary_disp_L_bp", "boundary_disp_R_bp", "S_mobility",
            "discovery_verdict", "longest_orf_bp", "orf_frac", "rrna_genes",
            "is_hmm_hit", "is_hmm_evalue", "domain_hit", "domain_evalue",
            "mobility_strong", "class_call"]
    out = open(os.path.join(args.outdir, "element_annotation.tsv"), "w")
    out.write("\t".join(cols) + "\n")

    tally = defaultdict(int)
    for nm in sorted(lens):
        L = lens[nm]
        st = struct.get(nm, {})
        o = orf.get(nm, 0)
        frac = o / L if L else 0.0
        rr = ("." if rrna is None else
              ",".join(sorted(set(rrna.get(nm, [])))) or ".")
        rr_unknown = rrna is None
        ishit, isev = best(ishits.get(nm) if ishits is not None else None)
        dhit, dev = best(dhits.get(nm) if dhits is not None else None)

        def fget(k, d=float("nan")):
            try:
                return float(st.get(k, ""))
            except (TypeError, ValueError):
                return d
        S = fget("S_mobility")
        ncop = fget("n_copies")
        strong = (S == S and S >= args.strong_s) or (ncop == ncop and ncop >= 5)
        verdict = st.get("verdict", ".")

        # --- derived class. Evidence columns above are the real output. ---
        if rr != "." and frac < args.min_orf_frac:
            cls = "ribosomal_repeat"
        elif ishit != ".":
            cls = "known_IS_missed_by_ISEScan"
        elif dhit != "." and IS_WORDS.search(dhit):
            cls = "known_IS_missed_by_ISEScan"
        elif dhit != "." and MGE_WORDS.search(dhit):
            cls = "known_other_MGE"
        elif verdict == "SEGMENTAL_DUP_LIKE":
            cls = "segmental_duplication"
        elif informative(dhit):
            # identified, just not as an MGE (Rhs toxin, LacI operon, ...)
            cls = "known_non_MGE_repeat"
        elif frac >= args.min_orf_frac and strong:
            cls = "unknown_mobile_element_candidate"
        elif strong:
            # ORF-poor but structurally mobile. This is a PROMOTION: the project
            # exists to find elements that are not conventional transposases.
            cls = "noncanonical_mobile_candidate"
        else:
            cls = "unclassified"
        if rr_unknown and cls == "ribosomal_repeat":
            cls = "unclassified"
        tally[cls] += 1

        out.write("\t".join(map(str, [
            nm, L, st.get("n_copies", "."), st.get("n_distinct_loci", "."),
            st.get("boundary_disp_L_bp", "."), st.get("boundary_disp_R_bp", "."),
            st.get("S_mobility", "."), verdict, o, "%.3f" % frac, rr,
            ishit, isev or ".", dhit, dev or ".", int(bool(strong)), cls])) + "\n")
    out.close()

    log("=" * 70)
    for k, v in sorted(tally.items(), key=lambda kv: -kv[1]):
        log("  %-34s %d" % (k, v))
    log("=" * 70)
    if rrna is None:
        log("NOTE: barrnap did not run -- the rRNA axis is UNKNOWN, not negative.")
    log("Flags only. No candidate was removed on family grounds.")
    log("wrote %s/element_annotation.tsv" % args.outdir)


if __name__ == "__main__":
    main()
