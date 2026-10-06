# FNA-only mobile-element discovery

Reference-free, annotation-free discovery of insertion elements by comparing
assembled genomes. No SRA download, no read alignment, no BAM.

**Scope: insertion elements under 5 kb.** Anything larger is out of the primary
catalogue by `armA.max_len` / `armB.max_insert` (both 5000) — but in Arm B it is
**parked, not discarded.** An IS110 nested inside a larger insertion is reported
as the larger event, and 2 of 13 IS110 events in the Arm B test were exactly
that (33 kb and 36 kb bubbles carrying an IS110 transposase). So oversized
bubbles go to `armB_large_events.tsv` + `armB_large_inserts.fna` with ids
`<cluster>.L%06d`, and `armB_events.tsv` carries a `parent_large_event_id`
column to link back. On the 3-genome test this retains 42 bubbles (median 15 kb,
range 5.5–66 kb) that the size gate used to throw away. Stage 32, which mines
them for the 1–5 kb repeated component inside the cargo, is not yet written.

> **Session continuity:** `docs/MEMORY_BACKUP.md` carries the working state
> (including the numbers below that this README now predates), the open
> items, and the operational gotchas. Re-sync it with
> `./tools/backup_memory.sh` after any session that changes memory.

Two arms, because they have complementary blind spots:

```
                         FNA collection
                               │
              ┌────────────────┴────────────────┐
              ▼                                 ▼
   ARM A  within-genome                ARM B  between-genome
   multi-copy repeats                  empty vs filled allele
              │                                 │
   copy1 AAAA─[═══]─CCCC              A: LEFT ─────────── RIGHT
   copy2 GGGT─[═══]─TTAT              B: LEFT ─ INSERT ── RIGHT
   copy3 CTAC─[═══]─AGGG              C: LEFT ─ INSERT ── RIGHT
         ↑ divergent flanks ↑         D: LEFT ─────────── RIGHT
              │                                 │
              └────────────────┬────────────────┘
                               ▼
                    universal element catalogue
                               │
                               ▼
                    ANNOTATION LAST  (stage 50)
              IS110 / IS30 / IS903 / prophage / unknown
```

**Annotation runs last, never first.** Nothing upstream of stage 50 knows what a
transposase is. That is the whole point: starting from an IS library can only
find what the library already contains.

---

## Pipeline

| stage | script | what it does | tool |
|---|---|---|---|
| 00 | (in stage 10) | assembly QC: N50, contigs, N fraction | — |
| 10 | `10_skani_cluster.py` | ANI clusters, so compared loci are homologous | **skani** |
| 20 | `20_armA_multicopy.py` | **Arm A**: multi-copy + boundary voting | minimap2 |
| 30 | `30_armB_graph.py` | **Arm B**: graph bubbles → empty/filled alleles | **minigraph + gfatools** |
| 30b | `30b_armB_pangraph.py` | Arm B, alternative engine | **PanGraph 1.4** |
| 31 | `31_armB_refine.py` | nucleotide junctions, insert sequence, **TSD** | **svim-asm** / nucdiff |
| 40 | `40_merge_catalog.py` | union both arms, dedup, presence/absence | mmseqs2 |
| 50 | `50_annotate.py` | rRNA screen, transposase HMM, IS family | barrnap, prodigal+HMMER, **ISEScan** |
| 60 | `60_benchmark_armA_vs_isescan.py` | is Arm A a real discovery engine? | ISEScan |
| 61 | `61_benchmark_vs_readgold.py` | FNA-only vs your read-validated junctions | — |

## Quick start

### What already works, with no conda solve

The graph tools are static binaries, already installed in `env/bin/`:

| binary | version |
|---|---|
| `minigraph` | 0.21-r606 |
| `gfatools` | built from source |
| `skani` | 0.3.2 |
| `pangraph` | 1.4.0 |

Combined with the existing `claude-env` (minimap2 2.30, MUMmer4, mmseqs2,
samtools, Biopython), **stages 10, 20, 30, 30b and 40 run today**:

```bash
export PATH=$PWD/env/bin:$HOME/.conda/envs/claude-env/bin:$PATH

# Arm A -- one genome, ~1.2 s, needs only minimap2 + stdlib
python3 scripts/20_armA_multicopy.py genome.fna --out out/g --threads 8

# Arm B -- ANI cluster, then graph bubbles
python3 scripts/10_skani_cluster.py --fna-list all_fna.txt --outdir out/stage10 --threads 16
python3 scripts/30_armB_graph.py --manifest out/stage10/clusters/cl0000.manifest \
        --outdir out/stage30/cl0000 --cluster-id cl0000 --threads 16

# union both arms
python3 scripts/40_merge_catalog.py \
   --armA-glob 'out/stage20/*_families.tsv' \
   --armA-elements-glob 'out/stage20/*_elements.fna' \
   --armB-events out/stage30/cl0000/armB_events.tsv \
   --armB-inserts out/stage30/cl0000/armB_inserts.fna --outdir out/stage40
```

Stage 30 emits the inserted sequence, the empty-site sequence and the TSD
straight from the graph alleles, so **stage 31 is an optional precision upgrade,
not a prerequisite** — stage 40 detects its absence and falls back to
graph-resolution coordinates.

### The remaining conda envs

`env/setup_env.sh` builds three envs, but a single mamba solve over ~18 packages
stalled on this filesystem for 25+ minutes at 0% CPU. Small solves are the fix,
and they only gate the two optional stages:

```bash
E=$HOME/.conda/envs
mamba create -y -p $E/fnains_sv    -c conda-forge -c bioconda svim-asm minimap2 samtools  # stage 31
mamba create -y -p $E/fnains_annot -c conda-forge -c bioconda isescan barrnap             # stage 50
mamba create -y -p $E/fnains_util  -c conda-forge -c bioconda seqkit nucdiff
```

`prodigal` and `hmmer` for stage 50 are already available as cluster modules
(`bio/prodigal/2.6.3`, `bio/hmmer/3.4`), so stage 50 can run its ORF and HMM
screens before ISEScan is installed.

### Full run

**There is no Snakemake entry point.** The Snakefile described the original
stage-10/40/50 pipeline, three of whose scripts were retired during the
simplification work and one of whose inputs (`stage40/mmseqs_rep_seq.fasta`)
no rule ever produced — it could not build a DAG on a fresh checkout. It is
kept as `scripts/_retired/Snakefile.legacy`; see the note beside it.

The pipeline that runs is the eleven stages in
`result_note/01_fna_to_insertion_call.md`. Every stage is submitted through
the generic wrapper, never from the login node:

```bash
sbatch --export=ALL,SCRIPT=scripts/100_event_dedup.py,ARGS="--recon 'l2/*/' --out dedup" \
       --output=$SCRATCH/logs/dedup_%j.out slurm/analysis.sh
```

Downloads are the one exception and run by hand on the transfer node
(`tools/download_on_transfer_node.sh`), never under sbatch.

Or on SLURM, per your existing conventions:

```bash
sbatch --export=FNA_LIST=/path/all_fna.txt,OUT=$SCRATCH/run01/armA slurm/armA_array.sh
sbatch --export=CLUSTERDIR=$SCRATCH/run01/stage10,OUT=$SCRATCH/run01/armB slurm/armB_array.sh
```

---

## What was cut from the original design, and why

The plan as written would have worked, but five pieces were paying no rent.

**1. Running both NucDiff *and* SVIM-asm — dropped one.**
They answer the same question (pairwise assembly-vs-assembly difference calling).
SVIM-asm is primary because it emits VCF with `POS`/`SVLEN`/inserted sequence, so
the existing Sniffles-shaped downstream flank code needs no changes. NucDiff is
kept only as the fragmented-draft fallback, where SVIM-asm's single-best-alignment
model is weakest.

**2. SyRI — demoted out of the default path.**
SyRI needs chromosome-level, largely one-to-one assemblies. Most of the staged
FNAs are short-read drafts, where it will either fail or return
synteny-block noise. It stays available for complete/long-read genomes only.

**3. RepeatModeler — dropped entirely.**
Built for eukaryotic TE families, costs hours per genome, and for bacterial
self-repeat discovery a direct minimap2 self-alignment does the same job in
about one second. Measured: 1.2 s wall for a complete 4.76 Mb *E. coli*.

**4. panISa — dropped.**
It is read/BAM-based, which is the exact dependency this pipeline exists to
remove.

**5. ANI clustering on identity alone — replaced.**
`--min-ani 99` is necessary but not sufficient. Two genomes can be 99% identical
over the 30% of themselves that aligns and share no architecture. The real guard
is **aligned fraction**, so `--min-af 85` is enforced alongside ANI.

### Added, because the design was missing them

**`minigraph` + `gfatools bubble` as the default Arm B engine.**
PanGraph is the better *concept* for this problem and is kept (v1.4.0, static
binary, no build). But its published scale is hundreds of genomes, and there are
13k staged FNAs. minigraph is built for exactly SV-bubble discovery, scales
further, and `gfatools bubble` reports the shortest and longest allele length per
bubble — the empty/filled length difference falls straight out. Both engines emit
the same event schema, so they are interchangeable and cross-checkable.

**Target-site duplication detection** (stages 30 and 31), as an *optional
annotation only*. Once empty and filled alleles are both in hand, TSD is free to
test:

```
empty  : LEFT [tsd] RIGHT
filled : LEFT [tsd] ELEMENT [tsd] RIGHT
```

**It is never a filter, and it must not be.** The IS110/IS1111 family leaves no
target-site duplication and has no terminal inverted repeats either. Measured on
this pipeline's own output: **0 of 13 IS110-positive events carry a TSD, against
43 of 170 for every other family.** Gating on TSD would discard 100% of IS110
events. So a TSD present is positive evidence of transposition; a TSD absent is
not evidence against it. Verified by grep that no caller in the pipeline branches
on `tsd_*`.

The same fact is why the *absence* of terminal inverted repeats does not hurt
this pipeline but does hurt ISEScan-style tools: Arm A keys on multi-copy
structure and flank divergence, neither of which requires a TIR.

**Boundary *voting* in a common frame** (stage 20). Merged self-alignment
intervals have ragged edges. Every copy is re-aligned to a representative copy,
so all boundary votes live in one coordinate frame; the consensus is the median
and its dispersion is reported in bp. On real *E. coli* IS this gives
**0.0–1.0 bp dispersion** — comparable to a read-based junction call.

---

## Validated behaviour

Everything below was actually run on staged *E. coli* assemblies, not projected.

### Arm A at scale -- 200 genomes, `--max-len 5000`, 0 failures

1.2 s per genome, no HMM, no IS library, no TIR requirement. 1,648 families:

| verdict | n | % |
|---|---:|---:|
| `STRONG_MOBILE_CANDIDATE` | 1,055 | 64.0% |
| `MOBILE_CANDIDATE` | 284 | 17.2% |
| `REPEAT_NO_BOUNDARY` (rrn-like) | 298 | 18.1% |
| `SEGMENTAL_DUP_LIKE` | 9 | 0.5% |
| `WEAK` | 2 | 0.1% |

STRONG candidates: median length **1,310 bp** (p10 768, p90 2,699), median
**4 distinct loci**, median internal identity **99.71%**.

**Boundary voting accuracy** (2,110 edges, cross-copy dispersion):

| ≤0 bp | ≤1 bp | ≤2 bp | ≤5 bp |
|---:|---:|---:|---:|
| 61.9% | 83.1% | **86.6%** | **89.3%** |

This is internal consistency between copies, *not* agreement with a read-based
gold standard — that is what `61_benchmark_vs_readgold.py` measures. But it shows
the FNA-only boundary is not the soft ±100 bp estimate one might expect from
merged alignment intervals.

The conservation profile behaves as a step function at the inferred edge
(family A0001, left boundary, cross-copy k-mer Jaccard):

```
offset  -500 -400 -300 -200 -100  -50 │   0  +50 +100 +150 +200
jaccard  .00  .00  .00  .00  .00  .00 │ .95  .95  .95  .95  .95
                       outside  ──────┼────── inside element
```

### IS110 as the worked example

IS110 is the hard case for conventional IS finders and the easy case for this
one. It has **no terminal inverted repeats** and leaves **no target-site
duplication**, so the two structural features most IS callers key on are both
absent. What it does have is high copy number per genome and sharp
flank divergence — exactly what Arm A keys on.

Scoring the 200-genome Arm A output against your own HMMs
(`hmm/is110_DEDD_Tnp20.hmm` = PF01548 DEDD + PF02371 Tnp20, both domains
required), with the HMM applied *only after* discovery:

| | IS110-positive families |
|---|---|
| families | **56**, in 55 of 200 genomes |
| length | median **1,310 bp** (range 1,293–3,989, all < 5 kb) |
| copies per genome | median **8** (range 3–23) — twice the overall STRONG median of 4 |
| internal identity | median 99.55% |
| boundary dispersion | median **0.0 bp** left, **1.0 bp** right |
| S_mobility | median **0.99** |
| verdict | 52 STRONG (93%), 3 MOBILE, **1 REPEAT_NO_BOUNDARY (missed)** |

So annotation-free discovery recovers IS110 at 93% STRONG with sub-2 bp boundary
agreement, and does it in 1.2 s per genome. The one `REPEAT_NO_BOUNDARY` miss and
one family with 859 bp left-boundary dispersion are nested/composite loci where
the merged interval spans more than the element.

**And the TSD result, which is the reason TSD is not a filter:**

| | events | with TSD |
|---|---:|---:|
| IS110-positive | 13 | **0 (0.0%)** |
| everything else | 170 | 43 (25.3%) |

IS110 insert sizes in Arm B were 1279 bp (×8) and 1448–1450 bp (×3), all with
zero TSD. Any pipeline that required a TSD would have thrown away every one of
them.

### Arm B -- 3-genome cluster, 248 bubbles → 183 empty/filled alleles

skani on 40 genomes gave 5 strain-level clusters (9/4/4/3/3) in 30 s. minigraph
graph build is the dominant cost: ~90 s for 3 *E. coli* genomes on 8 threads,
roughly linear in genome count. Genotyping is ~11 s per genome.

All 183 events came out `qc_tier: high`. Insert-size distribution again peaks in
the IS window (91 of 183 at 700–1,500 bp).

**TSD, where it exists, lands exactly where it should.** 43/183 events (23%)
carry a complexity-filtered TSD, and the length distribution is not flat:

```
TSD length:  3   4   5   6   7   8   9  10  12  13  15  28
     count:  2  13   4   5   5   1   7   1   2   1   1   1
                 ▲                   ▲
              IS3 family          IS1 / Tn5 / IS10
              (3-4 bp)            (9 bp)
```

Those two peaks are the canonical bacterial TSD lengths, and they were derived
without any IS model. The 777 bp inserts specifically carry 9 bp TSDs
(`TAATGATTT`, `CGCTTCCGG`, `GAATGAGTA`) — 777 bp ≈ IS1 at 768 bp, which makes a
9 bp TSD. Independent confirmation that those are real transposition events.

The other 77% are not failures. See the IS110 section below: the family this
project cares about most is in that 77% *by biology*, not by detection error.

### Stage 40 -- the union is doing real work

1,755 Arm A families + 183 Arm B events → **310 universal elements**:

| support | n |
|---|---:|
| `ARM_A_MULTICOPY` only | 202 |
| `ARM_B_EMPTY_FILLED` only | 87 |
| **`BOTH_ARMS`** | **21** |

87 elements were invisible to Arm A and 202 invisible to Arm B, which is the
argument for running both rather than picking one. Top `BOTH_ARMS` hits:

| element | len | max copies in one genome | genomes | S_mobility | TSD |
|---|---:|---:|---:|---:|---:|
| E00000 | 769 | 28 | 134 | 1.00 | **9 bp** |
| E00005 | 1329 | 27 | 56 | 1.00 | 10 bp |
| E00007 | 1310 | 30 | 48 | 1.00 | 7 bp |
| E00008 | 1345 | 9 | 50 | 1.00 | 12 bp |
| E00004 | 1259 | 7 | 74 | 0.97 | — |
| E00089 | 10976 | 3 | 2 | 0.50 | 9 bp |

E00000 — 769 bp, up to 28 copies per genome, present in 134/200 genomes, 9 bp
TSD — is IS1 by every structural criterion, identified with no reference to
IS1. E00089 at 11 kb with a 9 bp TSD is the large-element end of the same
signal.

### The honest false positive

Arm A family A0000 in the single-genome test: seven copies, 98.6% identity,
sharp edges after refinement — and it is the *rrn* operon. It carries 16S, 23S
and 5S motifs and its longest ORF is 252 bp, against 501–1,134 bp for the real
IS families in the same genome. Discovery cannot reject it without knowing what
rRNA is, so it must not try; stage 50 removes it with barrnap + ORF fraction,
which are structural screens rather than IS screens, so unbiasedness survives.
At scale this class is 15.8% of families and is already separated by
`REPEAT_NO_BOUNDARY` (median sharpness 0.00 vs 1.00 for STRONG).

## Known limits — stated, not papered over

**Polarity is undecidable from two assemblies.** `A = L-X-R` and `B = L-R` cannot
tell insertion from deletion. Stage 40 therefore emits `PRESENCE_ABSENCE`, and
only promotes to `GAIN`/`LOSS` when `--outgroup` supplies genomes outside the
cluster. This is an information limit, not a missing feature.

**Contig breaks cannot be repaired by algorithm.** Every copy and event carries a
`qc_tier`:

| tier | meaning |
|---|---|
| `high` | element + both ≥500 bp flanks inside one contig |
| `medium` | one side within 100–500 bp of a contig end |
| `unresolvable_contig_break` | the junction *is* the assembly gap |

Long-read/complete assemblies will substantially outperform short-read drafts
here, and no parameter changes that.

**TSD absence is not a negative result.** IS110/IS1111 leaves no TSD and no TIR
(0/13 measured). Nothing in the pipeline filters on `tsd_*`, and nothing should
be added that does. If you later add a scoring term for TSD, make it additive
positive evidence only — never a requirement.

**Arm A is not unbiased on its own.** It favours high copy number, recent
activity, and elements the assembler resolved; it misses singletons, diverged
copies, and collapsed repeats. That is precisely why Arm B runs alongside it and
stage 40 takes the union. `support` in `universal_elements.tsv` records which
arms saw each element — `BOTH_ARMS` is the strongest class: repeated *within* a
genome and polymorphic *between* genomes.

**ISCompare is a benchmark, never the discovery engine.** It finds IS first, then
compares location, so it carries the family bias this pipeline is built to avoid.
Cloned by `setup_env.sh` under `env/ext/` for positive-control use only.

---

## Key outputs

`universal_elements.tsv` — one row per element:
`element_id, len_median, n_armA_families, n_armB_events, support,
max_copies_in_one_genome, n_genomes_with_copy, presence_count, absence_count,
best_S_mobility, max_tsd_len, allele_state, gain_loss_confidence`

`universal_events.tsv` — one row per site. For Arm A each row is one insertion-site
context, which is the multi-site bag for downstream modelling: one element →
{site₁, site₂, site₃, …}, all within one genome, so GC background, species,
chromosome composition and assembly pipeline are held constant.

`S_mobility` (stage 20), components reported separately so it can be recalibrated:

```
S = 0.25·internal_conservation + 0.30·flank_divergence
  + 0.25·boundary_sharpness    + 0.20·locus_multiplicity
```

Verdicts are gated on boundary evidence, not on the score alone — a repeat with
no sharp edge becomes `REPEAT_NO_BOUNDARY` regardless of how high `S` climbs.

---

## Module status and capability boundary

Frozen 2026-08-28. Re-check with `scripts/73_layer2_regression.py` after any
Layer 2 change.

| layer | component | status |
|---|---|---|
| **1** | candidate discovery (Arm A / Arm B) | in use |
| **2** | `place_locus()` | **LOCKED** |
| **2** | exact decomposition | **LOCKED** |
| **2** | tolerant decomposition | **LOCKED** |
| **2** | placement-ambiguity reporting | **LOCKED** |
| **3** | polarity / ancestral-allele assignment | **NOT YET ASSIGNED** |
| — | nucleotide-junction certification | **NOT CLAIMED** |

**Layer 2 claims:** given a candidate locus and a cluster of homologous
assemblies, reconstruct the alternative alleles between homologous anchors,
preserve zero-length and retained target intervals, and decompose
high-homology long/short allele pairs into insertion/replacement
representations with explicit placement and junction ambiguity.

**Layer 2 does not claim** a nucleotide-resolved biochemical insertion
junction. Arm B's `junction_span` has median 35 bp (p75 282), so no
assembly-only comparison can certify ±2 bp; that requires read-level gold.
Validation is therefore two-tier:

| tier | gold | certifies |
|---|---|---|
| assembly | Arm B | detection, allele length, inserted sequence, coarse junction concordance |
| read | the 70 read-validated junctions | nucleotide junction, asymmetric local edit, microhomology |

The read tier has not been run — those junctions have never been located.

### Frozen Layer 2 benchmark (251 Arm B events, cl0000/cl0001/cl0002)

| gate | value | threshold |
|---|---:|---|
| resolved loci | 234 | ≥ 225 |
| catastrophic placements | 0 | = 0 |
| exact-method inserted length exact | 91.7% | ≥ 89% |
| tolerant-method inserted length exact | 74.0% | ≥ 71% |
| all decomposed | 82.0% | ≥ 79% |
| agreement regressions | 0 | = 0 |
| reverse-complement invariance | 96.2% | ≥ 95% |
| placement ambiguous | 4.3% | ≤ 8% |

Residual cross-method disagreements are recorded as
`validation_flag = assembly_crossmethod_discordant`, not resolved: neither
method is a length gold. Currently 1 in 234 loci (0.43%) with no repeated
failure signature. Revisit the algorithm only if that flag accumulates ~10+
cases sharing one shape.


### Layer 3 — outgroup architecture, CLOSED 2026-08-29

Re-check with `scripts/89_layer3_regression.py` after any Layer 3 change.

| panel | ANI to ingroup | role |
|---|---|---|
| **NEAR** | 98.2–99.0% | **primary inference** |
| FAR | 96.0–98.0% | perturbation validation only |
| MIXED | 3 near + 2 far | perturbation validation only |

**Taxa carry two independent qualifications.** `topology_informative` — its
flanks place, so it can resolve and root the local tree. `state_informative` —
it *also* matches a known allele confidently, so it can carry a Fitch character.
A distant outgroup is routinely the first without being the second; such taxa
enter the tree with `state=None`. A `novel_outgroup_allele` is topology-useful,
not useless. Before this separation such taxa were dropped from the tree
entirely, which halved the far panel's tree count and made a denominator
artefact look like a real monophyly improvement.

**Measured, cl0000, 94 loci, all panels run with role separation:**

| panel | P1 | P2 | trees built | mono-fail rate | no-state |
|---|---:|---:|---:|---:|---:|
| NEAR | 30 | 7 | 65 | 25% | 29 |
| FAR | 25 | 1 | 39 | 18% | 55 |
| MIXED | 31 | 3 | 60 | 27% | 34 |

MIXED gains one P1 and does **not** reduce monophyly failure. Intrusions fall
almost evenly on every panel member (near 16/16/15, far 15/15), so those loci
are non-monophyletic against essentially any outgroup — a property of the locus,
not of panel choice.

**The property that is frozen** is not coverage:

> Changing taxon sampling may make a P1 call disappear.
> It must never make a P1 call change its answer.

Across every panel tested, shared P1 loci agree on the ancestral allele —
19/19 and 25/25, **zero flips**.

**Persistent failures are kept, not rescued.** Loci non-monophyletic against any
panel are labelled `local_focal_clade_nonmonophyly` and stay unresolved, with
`monophyly_L/R/LR`, `outgroup_intruders_LR` and `mrca_support_LR` retained for
later analysis. They are deliberately **not** fixed by re-drawing the focal
clade: choosing a clade because it is monophyletic at this locus, then reading
ancestry off the same tree, is circular. The clade is defined once from
genome-wide ANI; the local genealogy tests whether it applies here.

**Not attempted:** `--allele-margin` stays frozen at 0.005. Distant outgroups
produce smaller margins, but that means the sequence no longer distinguishes the
alleles — not that the threshold should shrink. Loosening it is how a 55.7%
majority-vote flip rate would be reintroduced. Any future divergence-scaled
margin needs a truth-calibrated benchmark first.
