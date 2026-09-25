#!/usr/bin/env python3
"""Layer 3c -- polarity from LOCAL genealogy, not from counting outgroup votes.

STATUS: outgroup architecture CLOSED 2026-08-29. NEAR (98.2-99.0% ANI) is the
primary inference panel; FAR and MIXED panels are perturbation validation only
and must not be used to try to rescue monophyly. Re-verify any change with
scripts/89_layer3_regression.py.

TAXA HAVE TWO INDEPENDENT ROLES, and conflating them cost the far panel its
entire contribution until 2026-08-29:
    topology_informative -- flanks place, so it can resolve and root the tree
    state_informative    -- it also matches a known allele confidently
A distant outgroup is routinely the first without the second. Such taxa enter
the tree with state=None (fitch treats that as the universal set): useful for
topology, silent on the character. NOVEL is topology-useful, not useless.

The sensitivity grid killed majority voting: 34.4% of loci changed their
supported allele between max-ANI 99.0 and 98.5, and 7.8% flipped from a
UNANIMOUS call at 99.0 to a different real allele at 98.5. Close relatives can
all share the derived event, so "5/5 outgroups say B" is not ancestry.

What replaces it: for each locus, build the genealogy OF THAT LOCUS from its own
flanking sequence, then reconstruct the ingroup MRCA's allele on that tree.

    LEFT flank  -> MAFFT -> ML tree_L
    RIGHT flank -> MAFFT -> ML tree_R
    LEFT+RIGHT partitioned -> ML tree_LR   <- the tree polarity is read from

Why local flanks rather than whole-genome ANI: bacterial recombination means the
species tree need not be this locus's tree, and it is this locus's history that
decides whether the element was gained or lost. The variable interval itself is
EXCLUDED from the alignment -- including it would let the character being
reconstructed determine the tree that reconstructs it.

Why three trees: the left and right sides of a locus can have different
histories. Agreement between tree_L and tree_R is evidence the local genealogy
is real; disagreement is `local_genealogy_conflict`, and the locus is reported
unresolved rather than forced onto one concatenated tree.

Ancestral state: ML tree for the topology, Fitch parsimony for the discrete
allele character. Support is the fraction of UFBoot replicate trees giving the
same ingroup-MRCA state -- a bootstrap-backed number, not a vote count.

Hard gates, any of which forces unresolved:
  * ingroup not monophyletic at this locus (recombination, mis-placed locus,
    or a cluster definition that does not hold here);
  * tree_L and tree_R disagree on the reconstructed state;
  * fewer than --min-lineages informative independent outgroups;
  * Fitch returns an ambiguous state set at the ingroup MRCA.

Tiers:
  P1  all gates pass, bootstrap support >= --p1-support
  P2  supported but sampling-limited (one lineage, or moderate support)
  P0  unresolved

NJ is computed only as a debug sanity check and never enters the label.

    83_local_genealogy_polarity.py --loci loci.tsv --allele-table alleles.tsv \\
        --og-alleles og/cl0000/outgroup_alleles.tsv --seed-loci seed.tsv \\
        --seed-genomes cl0000.manifest --outgroups outgroups.manifest \\
        --outdir polarity/cl0000
"""
import argparse
import csv
import importlib.util
import os
import shutil
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib.fasta import Fasta, revcomp, write_fasta                # noqa: E402
from lib.util import log                                         # noqa: E402

_HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location(
    "_ar", os.path.join(_HERE, "70_allele_reconstructor.py"))
_ar = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_ar)

try:
    from Bio import Phylo
except ImportError:                                              # pragma: no cover
    Phylo = None


def sample_name(path):
    b = os.path.basename(path)
    for suf in (".gz", ".fna", ".fa", ".fasta"):
        if b.endswith(suf):
            b = b[: -len(suf)]
    return b


# ------------------------------------------------------------------ Fitch
def fitch(tree, states):
    """Fitch parsimony on a ROOTED Bio.Phylo tree.

    `states` maps leaf name -> state string, or None for unknown. Unknown
    leaves contribute the universal set, so they neither support nor block a
    reconstruction. Returns the state set at every clade.
    """
    alphabet = {s for s in states.values() if s}
    up = {}

    def post(cl):
        if cl.is_terminal():
            s = states.get(cl.name)
            up[cl] = {s} if s else set(alphabet)
            return up[cl]
        sets = [post(c) for c in cl.clades]
        inter = set.intersection(*sets) if sets else set()
        up[cl] = inter if inter else set().union(*sets) if sets else set()
        return up[cl]

    post(tree.root)
    return up


def mrca_state(tree, ingroup, states, outgroups=()):
    """Fitch state at the focal MRCA, monophyly, intruding taxa and support.

    `intruders` are taxa sitting INSIDE the focal clade's MRCA that do not
    belong to it. Naming them is what separates "this locus recombined" from
    "the clade still has an impurity": the same intruder on both sides points
    at the clade definition, different intruders point at recombination.
    """
    present = [n for n in ingroup if any(t.name == n for t in tree.get_terminals())]
    if len(present) < 2:
        return {"state": None, "mono": False, "intruders": (), "support": None}
    try:
        node = tree.common_ancestor(present)
    except (ValueError, LookupError):
        return {"state": None, "mono": False, "intruders": (), "support": None}
    under = {t.name for t in node.get_terminals()}
    intr = tuple(sorted(under - set(present)))
    # Not every intruder means the same thing. A NEAR-EXTERNAL genome -- a
    # discovery-cluster member that fell just below the all-pairs ANI floor --
    # sitting inside the focal MRCA is benign: it is ingroup-like, and its
    # presence only means the focal clade was drawn conservatively. What
    # actually invalidates polarity is a SELECTED OUTGROUP inside the clade,
    # because then the clade's ancestor is not ancestral to the outgroups and
    # the rooting that polarity depends on is gone.
    og_intr = tuple(x for x in intr if x in set(outgroups))
    up = fitch(tree, states)
    sup = None
    if node.confidence is not None:
        sup = node.confidence
    elif node.name:
        # IQ-TREE writes "SH-aLRT/UFboot" into the node label
        try:
            sup = float(str(node.name).split("/")[-1])
        except ValueError:
            sup = None
    return {"state": up.get(node), "mono": not og_intr, "strict_mono": not intr,
            "intruders": intr, "outgroup_intruders": og_intr, "support": sup}


def classify_monophyly(rL, rR, rLR, min_support=70.0):
    """Why did focal monophyly fail? Answer from topology, not from guessing."""
    mL, mR, mLR = rL["mono"], rR["mono"], rLR["mono"]
    if mL and mR and mLR:
        return "none"
    if mL and mR and not mLR:
        # both sides independently recover the clade; only the concatenation
        # does not. That is a concatenation artefact, NOT evidence the ingroup
        # is not a clade, and must not be lumped with real non-monophyly.
        return "concat_only"
    sup = [x for x in (rL["support"], rR["support"], rLR["support"])
           if x is not None]
    if sup and max(sup) < min_support:
        return "low_support"
    if not mL and mR:
        return "left_only"
    if not mR and mL:
        return "right_only"
    if not mL and not mR:
        return ("bilateral_same_lineage"
                if set(rL["intruders"]) & set(rR["intruders"])
                else "bilateral_different_lineages")
    return "other"


# ------------------------------------------------------------------ external
def run(cmd, **kw):
    return subprocess.run(cmd, stdout=subprocess.DEVNULL,
                          stderr=subprocess.DEVNULL, check=False, **kw).returncode


def mafft(inp, out, mafft_bin, threads=1):
    with open(out, "w") as fh:
        return subprocess.run([mafft_bin, "--auto", "--quiet",
                               "--thread", str(threads), inp],
                              stdout=fh, stderr=subprocess.DEVNULL,
                              check=False).returncode


def iqtree(aln, pre, iq, boot, threads=1, partition=None):
    cmd = [iq, "-s", aln, "-m", "MFP", "-T", str(threads), "--quiet", "-redo",
           "--prefix", pre]
    if boot:
        cmd += ["-B", str(boot), "--alrt", str(boot), "--wbtl"]
    if partition:
        cmd += ["-p", partition]
    return run(cmd)


def read_tree(path):
    if Phylo is None or not os.path.exists(path):
        return None
    try:
        return Phylo.read(path, "newick")
    except Exception:
        return None


def rooted(tree, outgroups):
    """Root with the outgroups so the ingroup MRCA has an ancestor."""
    if tree is None:
        return None
    names = {t.name for t in tree.get_terminals()}
    og = [o for o in outgroups if o in names]
    if not og:
        return tree
    try:
        tree.root_with_outgroup(*og) if len(og) > 1 else tree.root_with_outgroup(og[0])
    except Exception:
        try:
            tree.root_with_outgroup(og[0])
        except Exception:
            pass
    return tree


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--loci", required=True)
    ap.add_argument("--allele-table", required=True, help="Layer 2 alleles.tsv")
    ap.add_argument("--og-alleles", required=True, help="Layer 3b output")
    ap.add_argument("--seed-loci", required=True)
    ap.add_argument("--seed-genomes", required=True,
                    help="the full discovery cluster. Supplies anchor sequence "
                         "and contributes taxa to the tree. A genome only needs "
                         "to be here to define the locus; it does NOT have to "
                         "be in the focal clade.")
    ap.add_argument("--focal-clade", default=None,
                    help="polarity ingroup from 84_polarity_ingroup_refiner. "
                         "Monophyly and the MRCA are defined over THIS set. "
                         "Discovery-cluster genomes outside it stay in the tree "
                         "as near-external lineages -- they are real "
                         "observations and dropping them would discard "
                         "evidence -- but they do not define the ancestor. "
                         "Defaults to --seed-genomes.")
    ap.add_argument("--outgroups", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--flank", type=int, default=1000)
    ap.add_argument("--anchor", type=int, default=500)
    ap.add_argument("--max-interval", type=int, default=7000)
    ap.add_argument("--boot", type=int, default=1000)
    ap.add_argument("--min-lineages", type=int, default=2)
    ap.add_argument("--p1-support", type=float, default=0.80)
    ap.add_argument("--max-loci", type=int, default=0, help="0 = all")
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--iqtree", default="iqtree")
    ap.add_argument("--mafft", default="mafft")
    ap.add_argument("--minimap2", default="minimap2")
    ap.add_argument("--keep-tmp", action="store_true")
    args = ap.parse_args()

    if Phylo is None:
        log("FATAL: biopython Phylo unavailable")
        return 1
    os.makedirs(args.outdir, exist_ok=True)

    loci = {r["locus_id"]: r for r in
            csv.DictReader(open(args.loci), delimiter="\t")}
    seeds = {r["locus_id"]: r for r in
             csv.DictReader(open(args.seed_loci), delimiter="\t")}
    ing = {sample_name(p): p for p in
           (l.strip() for l in open(args.seed_genomes) if l.strip())}
    focal = set(ing)
    if args.focal_clade and os.path.exists(args.focal_clade):
        focal = {sample_name(p) for p in
                 (l.strip() for l in open(args.focal_clade) if l.strip())}
        focal &= set(ing)
    ogs = {sample_name(p): p for p in
           (l.strip() for l in open(args.outgroups) if l.strip())}
    log("focal clade: %d of %d discovery-cluster genomes"
        % (len(focal), len(ing)))

    # leaf allele states -------------------------------------------------
    ing_state = defaultdict(dict)
    for r in csv.DictReader(open(args.allele_table), delimiter="\t"):
        for g in r["genomes"].split(","):
            if g:
                ing_state[r["locus_id"]][g] = r["allele_id"]
    og_state = defaultdict(dict)
    for r in csv.DictReader(open(args.og_alleles), delimiter="\t"):
        if r["outgroup"] == ".":
            continue
        a = r["allele_assignment"]
        # only committed, non-novel assignments are informative characters
        if (r["allele_match_status"] == "unresolved" or a in (".", "", "NOVEL")
                or a.startswith("ambiguous_between")):
            continue
        og_state[r["locus_id"]][r["outgroup"]] = a

    lids = [l for l in sorted(loci) if l in seeds]
    if args.max_loci:
        lids = lids[: args.max_loci]
    log("%d loci, %d ingroup, %d outgroups" % (len(lids), len(ing), len(ogs)))

    fas = {}

    def fa(g, p):
        if g not in fas:
            fas[g] = Fasta(p)
        return fas[g]

    rows = []
    tiers = Counter()
    workroot = os.path.join(args.outdir, "trees")
    os.makedirs(workroot, exist_ok=True)

    for lid in lids:
        sd, L2 = seeds[lid], loci[lid]
        st = dict(ing_state.get(lid, {}))
        st.update(og_state.get(lid, {}))
        n_og_inf = len(og_state.get(lid, {}))
        base = {"locus_id": lid, "n_informative_outgroups": n_og_inf,
                "ingroup_monophyly": ".", "ancestral_allele": ".",
                "ancestral_bootstrap_support": "NA", "tree_L_state": ".",
                "tree_R_state": ".", "local_genealogy_conflict": ".",
                "polarity_status": "unresolved", "polarity_tier": "P0",
                "n_taxa": 0, "note": ""}
        if n_og_inf < args.min_lineages:
            base["note"] = "insufficient_informative_outgroups"
            rows.append(base)
            tiers["P0"] += 1
            continue

        # --- gather flanks, in each taxon's own placement of the locus ----
        tmp = tempfile.mkdtemp(prefix="gen.", dir=workroot)
        qp = os.path.join(tmp, "anchors.fna")
        gp0 = ing.get(sd["seed_genome"])
        if not gp0:
            shutil.rmtree(tmp, ignore_errors=True)
            base["note"] = "seed_genome_missing"
            rows.append(base)
            tiers["P0"] += 1
            continue
        f0 = fa(sd["seed_genome"], gp0)
        a, b = int(sd["start"]), int(sd["end"])
        if sd["contig"] not in f0 or a - args.anchor < 0 or \
                b + args.anchor > f0.length(sd["contig"]):
            shutil.rmtree(tmp, ignore_errors=True)
            base["note"] = "anchor_off_contig"
            rows.append(base)
            tiers["P0"] += 1
            continue
        with open(qp, "w") as qf:
            write_fasta(qf, "L", f0.fetch(sd["contig"], a - args.anchor, a))
            write_fasta(qf, "R", f0.fetch(sd["contig"], b, b + args.anchor))

        seqL, seqR = {}, {}
        # A taxon has TWO independent qualifications, and conflating them was
        # costing the far panel its whole contribution:
        #
        #   topology_informative -- its flanks place, so it can help resolve and
        #                           root the local tree;
        #   state_informative    -- it also matches a known ingroup allele
        #                           confidently, so it can carry a character.
        #
        # A distant outgroup is routinely the first without being the second: at
        # 96-98% ANI, 50.6% of placements are `ambiguous_between` and 15% are
        # `novel_outgroup_allele`. Excluding those from the tree threw away the
        # rooting evidence that made monophyly failures drop 16 -> 6. They now
        # enter the tree with state=None, which fitch() treats as the universal
        # set: uninformative for the character, but fully present for topology.
        # NOVEL is topology-useful and state-uninformative, not useless.
        for g, gp in list(ing.items()) + list(ogs.items()):
            paf = os.path.join(tmp, g + ".paf")
            # map_anchors returns the PAF path, not an exit code
            try:
                _ar.map_anchors(qp, gp, paf, 1, args.minimap2)
            except subprocess.CalledProcessError:
                continue
            hits = _ar.best_hits(paf, 90.0, 0.8)
            pl = _ar.place_locus(hits.get("L", []), hits.get("R", []),
                                 args.max_interval, None, 20000, 0.05)
            if pl is None or pl["status"] == "ambiguous":
                continue
            f = fa(g, gp)
            clen = f.length(pl["contig"])
            s, e = pl["start"], pl["end"]
            # flanks OUTSIDE the variable interval: the character being
            # reconstructed must not build the tree that reconstructs it
            lf = f.fetch(pl["contig"], max(0, s - args.flank), s)
            rf = f.fetch(pl["contig"], e, min(clen, e + args.flank))
            if pl["strand"] == "-":
                lf, rf = revcomp(rf), revcomp(lf)
            if len(lf) < args.flank * 0.6 or len(rf) < args.flank * 0.6:
                continue
            seqL[g], seqR[g] = lf, rf

        base["n_taxa"] = len(seqL)
        n_ing = sum(1 for g in seqL if g in focal)
        n_og = sum(1 for g in seqL if g in ogs)
        if n_ing < 2 or n_og < args.min_lineages:
            shutil.rmtree(tmp, ignore_errors=True)
            base["note"] = "too_few_taxa_with_flanks"
            rows.append(base)
            tiers["P0"] += 1
            continue

        # --- alignments and ML trees --------------------------------------
        ok = True
        for tag, seqs in (("L", seqL), ("R", seqR)):
            raw = os.path.join(tmp, tag + ".fna")
            with open(raw, "w") as fh:
                for g, s_ in sorted(seqs.items()):
                    write_fasta(fh, g, s_)
            if mafft(raw, os.path.join(tmp, tag + ".aln"), args.mafft) != 0:
                ok = False
        if not ok:
            shutil.rmtree(tmp, ignore_errors=True)
            base["note"] = "alignment_failed"
            rows.append(base)
            tiers["P0"] += 1
            continue

        # concatenated LEFT+RIGHT for the tree polarity is read from
        cat = os.path.join(tmp, "LR.aln")
        aln = {}
        for tag in ("L", "R"):
            cur = None
            for line in open(os.path.join(tmp, tag + ".aln")):
                if line.startswith(">"):
                    cur = line[1:].strip()
                    aln.setdefault(cur, {}).setdefault(tag, [])
                elif cur:
                    aln[cur][tag].append(line.strip())
        with open(cat, "w") as fh:
            for g in sorted(aln):
                if set(aln[g]) != {"L", "R"}:
                    continue
                write_fasta(fh, g, "".join(aln[g]["L"]) + "".join(aln[g]["R"]))

        states = {g: st.get(g) for g in seqL}      # None where uninformative
        n_topo = sum(1 for g in seqL if g in ogs)
        n_state = sum(1 for g in seqL if g in ogs and st.get(g))
        base["n_topology_informative_outgroups"] = n_topo
        base["n_state_informative_outgroups"] = n_state
        res = {}
        for tag, alnp in (("L", os.path.join(tmp, "L.aln")),
                          ("R", os.path.join(tmp, "R.aln")),
                          ("LR", cat)):
            pre = os.path.join(tmp, tag)
            iqtree(alnp, pre, args.iqtree, args.boot if tag == "LR" else 0,
                   args.threads)
            t = read_tree(pre + ".treefile")
            if t is None:
                res[tag] = {"state": None, "mono": False, "strict_mono": False,
                            "intruders": (), "outgroup_intruders": (),
                            "support": None}
                continue
            t = rooted(t, list(ogs))
            res[tag] = mrca_state(t, [g for g in seqL if g in focal], states,
                                  list(ogs))

        def one(x):
            s_ = x["state"]
            return sorted(s_)[0] if s_ and len(s_) == 1 else "."

        monoLR = res["LR"]["mono"]
        aLR, aL, aR = one(res["LR"]), one(res["L"]), one(res["R"])
        fclass = classify_monophyly(res["L"], res["R"], res["LR"])
        base["monophyly_L"] = "yes" if res["L"]["mono"] else "no"
        base["monophyly_R"] = "yes" if res["R"]["mono"] else "no"
        base["monophyly_LR"] = "yes" if monoLR else "no"
        base["monophyly_failure_class"] = fclass
        base["intruders_LR"] = ",".join(res["LR"]["intruders"]) or "."
        base["outgroup_intruders_LR"] = ",".join(
            res["LR"]["outgroup_intruders"]) or "."
        base["strict_monophyly_LR"] = "yes" if res["LR"]["strict_mono"] else "no"
        base["mrca_support_LR"] = ("%.1f" % res["LR"]["support"]
                                   if res["LR"]["support"] is not None else "NA")
        base["ingroup_monophyly"] = "yes" if monoLR else "no"
        base["tree_L_state"], base["tree_R_state"] = aL, aR
        conflict = (aL != "." and aR != "." and aL != aR)
        base["local_genealogy_conflict"] = "yes" if conflict else "no"
        base["ancestral_allele"] = aLR

        # --- bootstrap stability -----------------------------------------
        supp = float("nan")
        ub = os.path.join(tmp, "LR.ufboot")
        if os.path.exists(ub) and aLR != ".":
            hits = tot = 0
            try:
                for bt in Phylo.parse(ub, "newick"):
                    bt = rooted(bt, list(ogs))
                    rb = mrca_state(bt, [g for g in seqL if g in focal],
                                    states, list(ogs))
                    s_ = rb["state"]
                    if not s_:
                        continue
                    tot += 1
                    if len(s_) == 1 and sorted(s_)[0] == aLR:
                        hits += 1
            except Exception:
                tot = 0
            if tot:
                supp = hits / tot
        base["ancestral_bootstrap_support"] = ("%.3f" % supp) if supp == supp else "NA"

        # --- tier ----------------------------------------------------------
        if aLR == "." or not monoLR or conflict:
            tier, status = "P0", "unresolved"
            if not monoLR:
                # Named for what it is: at this locus the focal clade -- defined
                # independently from genome-wide ANI -- is not a clade. Measured
                # across NEAR/FAR/MIXED panels, intrusions fall almost evenly on
                # every panel member (near 16/16/15, far 15/15), so these loci
                # are non-monophyletic against essentially ANY outgroup. They are
                # a property of the locus, not of panel choice, and no outgroup
                # panel recovers them.
                #
                # They are deliberately NOT rescued by re-drawing the focal clade
                # to fit the local tree. Choosing a clade because it happens to be
                # monophyletic here, then reading ancestry off the same tree, is
                # circular. The clade is defined once, genome-wide; the local
                # genealogy tests whether it applies at this locus; when it does
                # not, the answer is unresolved.
                base["note"] = ("concat_only_monophyly_failure"
                                if fclass == "concat_only"
                                else "local_focal_clade_nonmonophyly")
            elif conflict:
                base["note"] = "left_right_genealogy_conflict"
            else:
                base["note"] = "ambiguous_ancestral_state"
        elif supp == supp and supp >= args.p1_support and \
                n_og_inf >= args.min_lineages:
            tier, status = "P1", "ancestral_supported"
        else:
            tier, status = "P2", "ancestral_candidate"
            base["note"] = "moderate_support_or_limited_sampling"
        base["polarity_tier"], base["polarity_status"] = tier, status
        tiers[tier] += 1
        rows.append(base)
        if not args.keep_tmp:
            shutil.rmtree(tmp, ignore_errors=True)

    cols = ["locus_id", "n_taxa", "n_informative_outgroups", "ingroup_monophyly",
            "monophyly_L", "monophyly_R", "monophyly_LR",
            "n_topology_informative_outgroups", "n_state_informative_outgroups",
            "monophyly_failure_class", "intruders_LR", "outgroup_intruders_LR",
            "strict_monophyly_LR", "mrca_support_LR",
            "tree_L_state", "tree_R_state", "local_genealogy_conflict",
            "ancestral_allele", "ancestral_bootstrap_support",
            "polarity_status", "polarity_tier", "note"]
    with open(os.path.join(args.outdir, "polarity.tsv"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")

    log("=" * 72)
    log("LOCAL-GENEALOGY POLARITY  (%d loci)" % len(rows))
    for t in ("P1", "P2", "P0"):
        log("  %-4s %5d  %5.1f%%" % (t, tiers[t],
                                     100.0 * tiers[t] / len(rows) if rows else 0))
    fc = Counter(r.get("monophyly_failure_class", "") for r in rows
                 if r.get("monophyly_failure_class") not in ("", "none", None))
    if fc:
        log("")
        log("  monophyly failure taxonomy:")
        for k, v in fc.most_common():
            log("    %-34s %4d" % (k, v))
    nt = Counter(r["note"] for r in rows if r["note"])
    if nt:
        log("")
        log("  reasons for non-P1:")
        for k, v in nt.most_common():
            log("    %-38s %4d" % (k, v))
    log("=" * 72)
    log("wrote %s/polarity.tsv" % args.outdir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
