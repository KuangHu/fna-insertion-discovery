#!/usr/bin/env python3
"""ARM B -- cross-genome empty-vs-filled allele discovery. Family-agnostic.

Builds a sequence-based pangenome graph over ONE ANI cluster, then reads the
bubbles. A bubble whose paths differ in length by >= --min-insert IS the
empty/filled structure:

    genome A:  LEFT ------------- RIGHT     short path  (empty allele)
    genome B:  LEFT -- INSERT --- RIGHT     long path   (filled allele)
    genome C:  LEFT -- INSERT --- RIGHT
    genome D:  LEFT ------------- RIGHT

Nothing here knows what an IS is. No transposase model, no IS library, no
annotation. Output is deliberately called a STRUCTURAL ALLELE, not an
insertion: with assemblies alone you cannot tell gain from loss without an
outgroup, so polarity is left to stage 40.

Engine: minigraph + gfatools (scales to hundreds of genomes per cluster).
See 30b_armB_pangraph.py for the PanGraph engine over the same schema.

Outputs (in --outdir):
  graph.gfa                the cluster graph
  bubbles.tsv              every bubble, with min/max allele length
  calls/<sample>.bed       per-sample allele at every bubble
  armB_events.tsv          bubbles that are empty/filled structural alleles
  armB_pairs.tsv           (filled, empty) genome pairs to refine in stage 31
"""
import argparse
import os
import statistics
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib.fasta import Fasta                                     # noqa: E402
from lib.tsd import diff_alleles, find_tsd                      # noqa: E402
from lib.util import log, run, seq_md5                                   # noqa: E402

MISSING = {".", "*", "", "?"}


def nonempty(path):
    """--reuse must not resurrect a truncated file from a killed run."""
    return os.path.exists(path) and os.path.getsize(path) > 0


def sample_name(path):
    b = os.path.basename(path)
    for suf in (".fna.gz", ".fa.gz", ".fasta.gz", ".fna", ".fa", ".fasta"):
        if b.endswith(suf):
            return b[:-len(suf)]
    return b


# ------------------------------------------------------------- graph building
def build_graph(manifest, gfa, threads, minigraph):
    """-cxggs = incremental graph construction; the FIRST fasta is the backbone."""
    paths = [l.strip() for l in open(manifest) if l.strip()]
    cmd = [minigraph, "-cxggs", "-t", str(threads)] + paths
    with open(gfa, "w") as fh:
        run(cmd, stdout=fh, capture_output=False)
    return paths


def bubbles(gfa, out_tsv, gfatools):
    """`gfatools bubble` emits 14 columns:

        0 chrom  1 start  2 end  3 n_seg  4 n_path  5 is_inversion
        6 len_shortest  7 len_longest  8-10 (-1 unless -d)  11 segment list
        12 shortest-path sequence      13 longest-path sequence

    Columns 6/7 are the empty and filled allele lengths, so their difference is
    the candidate insertion size -- and 12/13 hand over the two allele
    SEQUENCES directly, which makes stage 30 self-sufficient: the inserted
    sequence and its TSD can be read here, with stage 31 only sharpening the
    coordinates.
    """
    raw = out_tsv + ".raw"
    with open(raw, "w") as fh:
        run([gfatools, "bubble", gfa], stdout=fh, capture_output=False)
    rows = {}
    with open(raw) as fh, open(out_tsv, "w") as out:
        out.write("bubble_id\tchrom\tstart\tend\tn_seg\tn_path\tis_inversion\t"
                  "len_shortest\tlen_longest\tlen_delta\tsegments\n")
        for i, line in enumerate(fh):
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            try:
                chrom, st, en = f[0], int(f[1]), int(f[2])
                n_seg, n_path, inv = int(f[3]), int(f[4]), int(f[5])
                lo, hi = int(f[6]), int(f[7])
            except ValueError:
                continue
            bid = "b%06d" % i
            segs = f[11] if len(f) > 11 else "."
            s_short = f[12] if len(f) > 12 else ""
            s_long = f[13] if len(f) > 13 else ""
            rows[(chrom, st, en)] = (bid, inv, lo, hi, s_short, s_long)
            out.write("%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n"
                      % (bid, chrom, st, en, n_seg, n_path, inv, lo, hi, hi - lo, segs))
    return rows


# --------------------------------------------------------------- genotyping
def call_sample(gfa, fasta, out_bed, threads, minigraph):
    cmd = [minigraph, "-cxasm", "--call", "-t", str(threads), gfa, fasta]
    with open(out_bed, "w") as fh:
        run(cmd, stdout=fh, capture_output=False)
    return out_bed


def parse_call(bed):
    """`minigraph -cxasm --call` emits 6 columns:

        0 chrom  1 start  2 end  3 source node  4 sink node
        5 walk : allele_len : strand : query_contig : query_start : query_end

    e.g.  >s2:128:-:CP154596.1:4269163:4269298
    The allele length in field 5 is what genotypes empty vs filled; the query
    interval is what the contig-break QC needs. Parsing is positional but
    falls back to the trailing integers, so a format change degrades instead of
    crashing.
    Returns {(chrom, start, end): (allele_len, contig, qs, qe, strand, walk)}.
    """
    out = {}
    with open(bed) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            try:
                key = (f[0], int(f[1]), int(f[2]))
            except ValueError:
                continue
            rec = f[5].strip()
            if rec in MISSING:
                out[key] = (None, None, None, None, ".", ".")
                continue
            tok = rec.split(":")
            walk, alen, contig, qs, qe, strand = tok[0], None, None, None, None, "."
            if len(tok) >= 6:
                try:
                    alen = int(tok[1])
                    strand = tok[2]
                    contig = tok[3]
                    qs, qe = int(tok[4]), int(tok[5])
                except ValueError:
                    alen = None
            if alen is None:
                ints = [int(t) for t in tok if t.lstrip("-").isdigit()]
                alen = abs(ints[0]) if ints else None
            out[key] = (alen, contig, qs, qe, strand, walk)
    return out


# ------------------------------------------------------------- classification
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", required=True,
                    help="cluster manifest from stage 10 (backbone on line 1)")
    ap.add_argument("--min-genome-frac", type=float, default=0.5,
                    help="reject a manifest entry smaller than this fraction of "
                         "the panel's MEDIAN genome size. Relative, not "
                         "absolute: the old 1 MB floor was calibrated on 2-5 Mb "
                         "genomes and rejects whole species outright -- "
                         "M. pneumoniae is 0.83 Mb, C. trachomatis 1.07, "
                         "T. pallidum 1.15. A median-relative floor catches what "
                         "the absolute one was for (plasmid-only GenBank entries "
                         "flagged 'Complete Genome': 55 in S. enterica, 7 in "
                         "P. aeruginosa) in every species regardless of size.")
    ap.add_argument("--min-genome-bytes", type=int, default=100000,
                    help="absolute floor as well, for a file that is not an "
                         "assembly at all (a 2,015-byte insert sequence once "
                         "passed as a genome and produced a silent empty run)")
    ap.add_argument("--allow-zero-bubbles", action="store_true",
                    help="record a genuine zero-bubble panel instead of "
                         "failing. Default is to FAIL: a silently failed panel "
                         "and a genuinely invariant one look identical.")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--cluster-id", default=None)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--min-insert", type=int, default=500,
                    help="min allele-length difference to call a structural allele")
    ap.add_argument("--max-insert", type=int, default=5000,
                    help="insertion elements only; >5 kb is out of scope")
    ap.add_argument("--len-tol", type=int, default=50,
                    help="allele lengths within this are the same allele")
    ap.add_argument("--min-empty", type=int, default=1)
    ap.add_argument("--min-filled", type=int, default=1)
    ap.add_argument("--flank-qc", type=int, default=500,
                    help="required unique anchor each side, inside one contig")
    ap.add_argument("--clean-span-tol", type=int, default=20,
                    help="empty-frame junction spans within this are at the "
                         "graph's resolution limit and are not called asymmetric")
    ap.add_argument("--min-tsd", type=int, default=3)
    ap.add_argument("--max-tsd", type=int, default=30)
    ap.add_argument("--keep-inversions", action="store_true")
    ap.add_argument("--minigraph", default="minigraph")
    ap.add_argument("--gfatools", default="gfatools")
    ap.add_argument("--reuse", action="store_true", help="reuse graph/calls if present")
    args = ap.parse_args()

    cid = args.cluster_id or os.path.basename(args.manifest).split(".")[0]
    os.makedirs(os.path.join(args.outdir, "calls"), exist_ok=True)
    gfa = os.path.join(args.outdir, "graph.gfa")

    paths = [l.strip() for l in open(args.manifest) if l.strip()]

    # ---- ASSERTION 1: every manifest entry must be a real assembly ---------
    # 2026-09-02. A test run listed 9 genomes of which 8 did not exist on disk
    # and the 9th was a 2,015-byte file containing an INSERT sequence
    # (">cl0001.B000082 GGAAGG...") that had been mis-saved under a genome
    # name -- an output fed back in as an input. minigraph got nothing, no
    # graph was written, all 9 call BEDs came out 0 bytes, and the stage still
    # wrote header-only TSVs and exited 0. At 234 loci a human noticed. Across
    # hundreds of programmatically generated panels, a silently failed panel
    # and a panel with genuinely no insertions are INDISTINGUISHABLE in the
    # output. Silent negative results are a forbidden failure mode.
    missing = [p for p in paths if not os.path.exists(p)]
    if missing:
        log("FATAL: %d of %d manifest entries do not exist, first: %s"
            % (len(missing), len(paths), missing[0]))
        return 2
    sizes = [os.path.getsize(p) for p in paths]
    med = sorted(sizes)[len(sizes) // 2]
    floor = max(args.min_genome_bytes, int(args.min_genome_frac * med))
    tiny = [(p, z) for p, z in zip(paths, sizes) if z < floor]
    if tiny:
        log("FATAL: %d of %d manifest entries are below %d bytes (%.0f%% of the "
            "panel median %d). These are not whole assemblies -- most often a "
            "plasmid-only GenBank entry flagged 'Complete Genome'. First: %s "
            "(%d bytes)" % (len(tiny), len(paths), floor,
                            100 * args.min_genome_frac, med, tiny[0][0], tiny[0][1]))
        return 2
    log("size check: median %.2f Mb, floor %.2f Mb, all %d entries pass"
        % (med / 1e6, floor / 1e6, len(paths)))
    log("%s: %d genomes, backbone=%s" % (cid, len(paths), sample_name(paths[0])))
    # bubble coordinates are in BACKBONE space, so the backbone supplies the
    # true flanking sequence -- the bubble's internal allele sequences do not
    # reach far enough out to see a TSD sitting against the source node.
    backbone = Fasta(paths[0])

    if not (args.reuse and nonempty(gfa)):
        log("building graph ...")
        build_graph(args.manifest, gfa, args.threads, args.minigraph)
    # ---- ASSERTION 2: the graph must actually have been built -------------
    if not nonempty(gfa):
        log("FATAL: %s is missing or empty after graph construction" % gfa)
        return 2
    n_seg = 0
    with open(gfa) as fh:
        for line in fh:
            if line.startswith("S"):
                n_seg += 1
                if n_seg >= 1:
                    break
    if n_seg == 0:
        log("FATAL: %s contains no S lines -- minigraph produced no graph" % gfa)
        return 2
    log("graph: %.1f MB" % (os.path.getsize(gfa) / 1e6))

    btsv = os.path.join(args.outdir, "bubbles.tsv")
    bindex = bubbles(gfa, btsv, args.gfatools)
    log("%d bubbles" % len(bindex))

    # genotype every genome against the graph
    calls, fais = {}, {}
    for p in paths:
        sn = sample_name(p)
        bed = os.path.join(args.outdir, "calls", sn + ".bed")
        if not (args.reuse and nonempty(bed)):
            call_sample(gfa, p, bed, args.threads, args.minigraph)
        calls[sn] = parse_call(bed)
        try:
            fais[sn] = Fasta(p)
        except Exception:
            fais[sn] = None
        log("  called %s (%d bubbles covered)"
            % (sn, sum(1 for v in calls[sn].values() if v[0] is not None)))

    ev = open(os.path.join(args.outdir, "armB_events.tsv"), "w")
    pr = open(os.path.join(args.outdir, "armB_pairs.tsv"), "w")
    al = open(os.path.join(args.outdir, "armB_inserts.fna"), "w")
    es = open(os.path.join(args.outdir, "armB_empty_sites.fna"), "w")
    # Oversized bubbles are OUT OF THE PRIMARY CATALOGUE, NOT DISCARDED. A large
    # cargo (prophage, ICE, plasmid segment) routinely carries a 1-5 kb element
    # inside it: 2 of 13 IS110 Arm B events were nested in 33 kb and 36 kb
    # bubbles. Dropping the bubble at the size gate would systematically lose
    # exactly the nested cases, so they are parked here for stage 32 to mine.
    lg = open(os.path.join(args.outdir, "armB_large_events.tsv"), "w")
    li = open(os.path.join(args.outdir, "armB_large_inserts.fna"), "w")
    lg.write("\t".join([
        "large_event_id", "cluster_id", "chrom", "bubble_start", "bubble_end",
        "is_inversion", "len_empty", "len_filled", "insert_size",
        "n_empty", "n_filled", "n_missing", "insert_md5",
        "empty_genomes", "filled_genomes"]) + "\n")
    ev.write("\t".join([
        "event_id", "cluster_id", "chrom", "bubble_start", "bubble_end",
        "is_inversion", "len_empty", "len_filled", "insert_size",
        "n_empty", "n_filled", "n_missing", "allele_state",
        "empty_ref_genome", "empty_ref_contig", "empty_ref_strand",
        "junction_L", "junction_R", "junction_span", "architecture",
        "insert_len_direct", "insert_md5", "tsd_len", "tsd_seq", "tsd_side",
        "tsd_conf", "empty_genomes", "filled_genomes", "qc_tier",
        "parent_large_event_id", "flags"]) + "\n")
    pr.write("event_id\tcluster_id\tfilled_genome\tempty_genome\tchrom\t"
             "bubble_start\tbubble_end\tinsert_size\n")

    def contiguity(sn):
        fa = fais.get(sn)
        return -(fa.total() / max(1, len(fa.names()))) if fa else 0

    n_ev = n_large = 0
    for key, (bid, inv, lo, hi, s_short, s_long) in bindex.items():
        chrom, bs, be = key
        if inv and not args.keep_inversions:
            continue
        obs = {}
        n_missing = 0
        for sn in calls:
            alen, contig, qs, qe, strand, _walk = calls[sn].get(key, (None,) * 6)
            if alen is None:
                n_missing += 1
            else:
                obs[sn] = (alen, contig, qs, qe, strand)
        if len(obs) < 2:
            continue
        lens = sorted(v[0] for v in obs.values())
        shortest = lens[0]
        empty = [s for s, v in obs.items() if v[0] <= shortest + args.len_tol]
        filled = [s for s, v in obs.items()
                  if v[0] >= shortest + args.min_insert]
        if len(empty) < args.min_empty or len(filled) < args.min_filled:
            continue
        len_filled = statistics.median([obs[s][0] for s in filled])
        insert = int(len_filled - statistics.median([obs[s][0] for s in empty]))
        if insert < args.min_insert:
            continue
        if insert > args.max_insert:
            # Out of the primary catalogue, but kept: stage 32 searches these
            # for the 1-5 kb repeated component sitting inside the large cargo.
            lseq, _, _ = diff_alleles(s_short, s_long)
            lid = "%s.L%06d" % (cid, n_large)
            lg.write("\t".join(map(str, [
                lid, cid, chrom, bs, be, inv,
                shortest, int(len_filled), insert,
                len(empty), len(filled), n_missing,
                seq_md5(lseq) if lseq else "NA",
                ",".join(sorted(empty)[:20]), ",".join(sorted(filled)[:20]),
            ])) + "\n")
            if lseq:
                li.write(">%s len=%d\n%s\n" % (lid, len(lseq), lseq))
            n_large += 1
            continue

        # FNA-specific QC: does a single contig actually span the locus?
        flags, tiers = [], []
        for s in filled + empty:
            alen, contig, qs, qe = obs[s][:4]
            fa = fais.get(s)
            if fa is None or contig is None or contig not in fa:
                tiers.append("medium")
                continue
            clen = fa.length(contig)
            left, right = min(qs, qe), clen - max(qs, qe)
            if left >= args.flank_qc and right >= args.flank_qc:
                tiers.append("high")
            elif left >= 100 and right >= 100:
                tiers.append("medium")
            else:
                tiers.append("unresolvable_contig_break")
        tier = ("high" if all(t == "high" for t in tiers)
                else "unresolvable_contig_break"
                if any(t == "unresolvable_contig_break" for t in tiers)
                else "medium")
        if n_missing:
            flags.append("missing_in_%d" % n_missing)
        if inv:
            flags.append("inversion")

        eid = "%s.B%06d" % (cid, n_ev)
        # the two allele sequences come straight out of gfatools bubble, so the
        # inserted sequence and its TSD are available without any realignment.
        # TSD is annotation, not a gate: IS110/IS1111 leaves none (0/13 measured
        # here, vs 43/170 for other families), so absence means nothing.
        ins_seq, aL, aR = diff_alleles(s_short, s_long)
        lflank = rflank = ""
        if chrom in backbone:
            clen_b = backbone.length(chrom)
            lflank = backbone.fetch(chrom, max(0, bs - args.flank_qc), bs)
            rflank = backbone.fetch(chrom, be, min(clen_b, be + args.flank_qc))
        tl, ts, tside, tconf = find_tsd(ins_seq,
                                       lflank + s_short[:aL],
                                       s_short[aL:] + rflank,
                                       min_tsd=args.min_tsd, max_tsd=args.max_tsd)
        # The two junctions are separate quantities and must stay separate.
        # diff_alleles anchored the insert inside the bubble, so:
        #   junction_L = where the filled allele departs from the empty one
        #   junction_R = where it rejoins
        # span == 0  -> clean insertion (both junctions coincide)
        # span  > 0  -> target-site deletion / asymmetric junction: the filled
        #               allele replaced `span` bp of the empty allele
        # Collapsing these to one POS is exactly what loses asymmetric events.
        # The junction pair is only meaningful in the frame of a genome that
        # actually LACKS the element. Bubble coordinates are backbone
        # coordinates, and the backbone is often a filled carrier -- reading
        # junctions off it then returns the element's own extent (spans of 777 /
        # 1349 bp, i.e. element-sized, which is the tell-tale of that mistake).
        # --call gives every genome its own interval at this bubble, so the
        # designated empty carrier supplies the correct frame.
        eref = min(empty, key=contiguity)
        _, e_contig, e_qs, e_qe, e_strand = obs[eref]
        if e_qs is None:
            jL = jR = -1
            span, arch = -1, "no_empty_frame"
        else:
            lo_q, hi_q = min(e_qs, e_qe), max(e_qs, e_qe)
            if e_strand == "-":
                # walk runs opposite the contig: the graph-left anchor trims the
                # contig-right end
                jL, jR = lo_q + aR, hi_q - aL
            else:
                jL, jR = lo_q + aL, hi_q - aR
            span = jR - jL
            # A span of a few bp at graph resolution cannot be told apart from
            # segment-boundary slack: minigraph's source/sink nodes do not break
            # exactly at the insertion point. So anything within tolerance is
            # reported as provisional, NOT as a confirmed asymmetric junction.
            # Only stage 31's pairwise nucleotide alignment can settle that.
            # Keep this field CATEGORICAL. The magnitude already has a column
            # of its own (junction_span); baking it into the label made every
            # asymmetric event its own singleton class and broke grouping in
            # the stage 61 architecture test.
            if abs(span) <= args.clean_span_tol:
                arch = "clean_or_graph_slack"
            else:
                arch = "asymmetric_candidate"
        ev.write("\t".join(map(str, [
            eid, cid, chrom, bs, be, inv,
            shortest, int(len_filled), insert,
            len(empty), len(filled), n_missing,
            "PRESENCE_ABSENCE",            # polarity decided in stage 40
            eref, e_contig or ".", e_strand or ".", jL, jR, span, arch,
            len(ins_seq), seq_md5(ins_seq) if ins_seq else "NA",
            tl, ts or ".", tside, tconf,
            ",".join(sorted(empty)[:20]), ",".join(sorted(filled)[:20]),
            tier, ".", ";".join(flags) or ".",
        ])) + "\n")
        if ins_seq:
            # header is exactly the event_id: stage 40 joins the catalogue on it
            al.write(">%s len=%d tsd=%d\n%s\n" % (eid, len(ins_seq), tl, ins_seq))
        es.write(">%s len=%d\n%s\n"
                 % (eid, len(lflank) + len(s_short) + len(rflank),
                    lflank + s_short + rflank))

        # refinement pairs: most contiguous filled x most contiguous empty
        for f in sorted(filled, key=contiguity)[:args.min_filled or 1]:
            for e in sorted(empty, key=contiguity)[:args.min_empty or 1]:
                pr.write("%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\n"
                         % (eid, cid, f, e, chrom, bs, be, insert))
        n_ev += 1

    ev.close()
    pr.close()
    al.close()
    es.close()
    lg.close()
    li.close()
    log("%s: %d empty/filled structural alleles" % (cid, n_ev))
    log("%s: %d oversized bubbles parked in armB_large_events.tsv (>%d bp, "
        "kept for nested sub-5kb search, NOT discarded)"
        % (cid, n_large, args.max_insert))

    # ---- ASSERTION 3: zero bubbles is a failure until proven otherwise -----
    # A panel with genuinely no length-polymorphic bubbles is possible in
    # principle -- a set of near-clonal assemblies can have none -- so this is
    # switchable. But the DEFAULT must be failure, because a silently failed
    # panel and a genuinely empty one are identical in the output, and at
    # census scale nobody inspects them one by one.
    if len(bindex) == 0 and not args.allow_zero_bubbles:
        log("FATAL: %s produced 0 bubbles from %d genomes. Either the panel is "
            "genuinely invariant or the run failed silently -- these look the "
            "same in the output. Re-run with --allow-zero-bubbles to record a "
            "genuine zero." % (cid, len(paths)))
        return 3
    return 0


if __name__ == "__main__":
    sys.exit(main())
