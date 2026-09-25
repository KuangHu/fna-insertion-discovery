# FNA-only mobile-element discovery.
#
#   snakemake -s Snakefile --configfile config/config.yaml -j 40 \
#       --config fna_list=/path/all_fna.txt run=run01
#
# Arm A is per-genome and ~1 s each, so it is left as one rule over the whole
# list rather than 13k jobs. Arm B fans out per ANI cluster.
import os

WORK = config["workdir"]
RUN  = config.get("run", "run01")
OUT  = os.path.join(WORK, RUN)
REPO = os.path.dirname(os.path.abspath(workflow.snakefile))
CORE  = config["envs"]["core"]
SVENV = config["envs"]["sv"]
ANENV = config["envs"]["annot"]
FNA_LIST = config["fna_list"]

def cids():
    d = os.path.join(OUT, "stage10", "clusters")
    if not os.path.isdir(d):
        return []
    return sorted(f[:-9] for f in os.listdir(d) if f.endswith(".manifest"))

rule all:
    input:
        os.path.join(OUT, "stage50", "element_annotation.tsv")

# ------------------------------------------------------------------ stage 10
rule cluster:
    input:  FNA_LIST
    output: tsv = os.path.join(OUT, "stage10", "clusters.tsv")
    threads: 16
    params: c = config["cluster"], d = os.path.join(OUT, "stage10")
    shell:
        "PATH={CORE}/bin:$PATH python3 {REPO}/scripts/10_skani_cluster.py "
        "--fna-list {input} --outdir {params.d} --threads {threads} "
        "--min-ani {params.c[min_ani]} --min-af {params.c[min_af]} "
        "--min-cluster-size {params.c[min_cluster_size]} "
        "--max-cluster-size {params.c[max_cluster_size]} "
        "--min-n50 {params.c[min_n50]} --max-contigs {params.c[max_contigs]}"

# ------------------------------------------------- stage 20 (Arm A per genome)
checkpoint armA:
    input:  FNA_LIST
    output: done = os.path.join(OUT, "stage20", "ARMA_DONE")
    threads: 40
    params: a = config["armA"], d = os.path.join(OUT, "stage20")
    shell:
        r"""
        mkdir -p {params.d}
        export PATH={CORE}/bin:$PATH
        cat {input} | xargs -P {threads} -I@ bash -c '
          f=@; s=$(basename "$f"); s=${{s%.gz}}; s=${{s%.fna}}; s=${{s%.fa}}; s=${{s%.fasta}}
          [ -s "{params.d}/$s"_families.tsv ] && exit 0
          w="$f"; case "$f" in *.gz) w="{params.d}/$s.fna"; zcat "$f" > "$w";; esac
          python3 {REPO}/scripts/20_armA_multicopy.py "$w" --out "{params.d}/$s" \
            --sample-id "$s" --threads 1 \
            --min-identity {params.a[min_identity]} --min-len {params.a[min_len]} \
            --max-len {params.a[max_len]} --min-copies {params.a[min_copies]} \
            --min-locus-sep {params.a[min_locus_sep]} --flank {params.a[flank]} \
            --kmer {params.a[kmer]} --min-sharpness {params.a[min_sharpness]} \
            --max-flank-homology {params.a[max_flank_homology]} \
            >> "{params.d}/$s.log" 2>&1 || echo "FAIL $s" >> {params.d}/failures.txt
          case "$f" in *.gz) rm -f "$w" "$w".fai;; esac'
        touch {output.done}
        """

# ------------------------------------------------- stage 30/31 (Arm B per cluster)
rule armB_cluster:
    input:  clusters = os.path.join(OUT, "stage10", "clusters.tsv")
    output: ev = os.path.join(OUT, "stage30", "{cid}", "armB_events.tsv"),
            rf = os.path.join(OUT, "stage30", "{cid}", "armB_refined.tsv")
    threads: 16
    params: b = config["armB"],
            man = lambda w: os.path.join(OUT, "stage10", "clusters", w.cid + ".manifest"),
            d = lambda w: os.path.join(OUT, "stage30", w.cid)
    shell:
        "PATH={CORE}/bin:$PATH python3 {REPO}/scripts/30_armB_graph.py "
        "--manifest {params.man} --outdir {params.d} --cluster-id {wildcards.cid} "
        "--threads {threads} --reuse "
        "--min-insert {params.b[min_insert]} --max-insert {params.b[max_insert]} "
        "--len-tol {params.b[len_tol]} --flank-qc {params.b[flank_qc]} && "
        "PATH={SVENV}/bin:$PATH python3 {REPO}/scripts/31_armB_refine.py "
        "--pairs {params.d}/armB_pairs.tsv --clusters {input.clusters} "
        "--outdir {params.d} --threads {threads} "
        "--min-tsd {params.b[min_tsd]} --max-tsd {params.b[max_tsd]}"

def armB_all(wildcards):
    checkpoints.armA.get(**wildcards)
    return expand(os.path.join(OUT, "stage30", "{cid}", "armB_refined.tsv"),
                  cid=cids())

# ------------------------------------------------------------------ stage 40
rule merge:
    input:  armB = armB_all,
            done = os.path.join(OUT, "stage20", "ARMA_DONE")
    output: cat = os.path.join(OUT, "stage40", "universal_elements.tsv"),
            fna = os.path.join(OUT, "stage40", "all_elements.fna")
    threads: 8
    params: m = config["merge"], d = os.path.join(OUT, "stage40"),
            og = lambda w: ("--outgroup " + config["merge"]["outgroup"]
                            if config["merge"].get("outgroup") else "")
    shell:
        r"""
        PATH={CORE}/bin:$PATH python3 {REPO}/scripts/40_merge_catalog.py \
          --armA-glob '{OUT}/stage20/*_families.tsv' \
          --armA-elements-glob '{OUT}/stage20/*_elements.fna' \
          --armB-refined <(head -1 $(ls {OUT}/stage30/*/armB_refined.tsv | head -1); \
                           tail -qn +2 {OUT}/stage30/*/armB_refined.tsv) \
          --armB-events  <(head -1 $(ls {OUT}/stage30/*/armB_events.tsv | head -1); \
                           tail -qn +2 {OUT}/stage30/*/armB_events.tsv) \
          --armB-inserts <(cat {OUT}/stage30/*/armB_inserts.fna) \
          --outdir {params.d} --threads {threads} \
          --min-seq-id {params.m[min_seq_id]} --cov {params.m[cov]} {params.og}
        """

# ---------------------------------------------- stage 50 (annotation LAST)
rule annotate:
    input:  cat = os.path.join(OUT, "stage40", "universal_elements.tsv"),
            fna = os.path.join(OUT, "stage40", "mmseqs_rep_seq.fasta")
    output: os.path.join(OUT, "stage50", "element_annotation.tsv")
    threads: 8
    params: a = config["annotate"], d = os.path.join(OUT, "stage50"),
            hmm = lambda w: ("--hmm " + config["annotate"]["hmm"]
                             if config["annotate"].get("hmm") else "")
    shell:
        "PATH={ANENV}/bin:$PATH python3 {REPO}/scripts/50_annotate.py "
        "--elements {input.fna} --catalog {input.cat} --outdir {params.d} "
        "--threads {threads} --min-orf-frac {params.a[min_orf_frac]} {params.hmm}"
