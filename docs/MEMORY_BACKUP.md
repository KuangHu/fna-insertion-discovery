# Memory backup — session continuity notes

The Claude Code memory for this project lives outside the repo, at:

```
~/.claude/projects/-global-home-users-kh36969-fna-based-mgefinder-project/memory/
```

That directory is **not** backed up and **not** in version control, and it has
already been lost once: the 2026-08-27 session that built and validated most of
this pipeline vanished with no transcript and no memory files, and the state had
to be reconstructed by reading file mtimes. This file is the durable copy. It
lives in the repo so it survives anything that clears `~/.claude`.

## Restore

```bash
M=~/.claude/projects/-global-home-users-kh36969-fna-based-mgefinder-project/memory
mkdir -p "$M" && cp docs/memory/*.md "$M"/
```

## Re-sync after changing memory

```bash
./tools/backup_memory.sh          # copies memory/ -> docs/memory/
```

Run it at the end of any session that added or edited a memory file.

---

## Snapshot of the state, 2026-08-27

Kept here in prose as well as in `docs/memory/`, so this file is useful even to
a reader who never restores the memory directory.

### What the pipeline is

Reference-free insertion-element discovery from assemblies only — no SRA, no
read alignment, no BAM. Two arms with complementary blind spots, unioned:

- **Arm A** — within-genome multi-copy repeats with divergent flanks.
- **Arm B** — between-genome empty vs filled alleles from a minigraph bubble.

Annotation runs **last** (stage 50). Primary catalogue is 300 bp – 5 kb.
`README.md` holds the full design rationale and is reliable **except** for the
Arm B and stage 40 numbers, which it predates — see the table below.

**This architecture is settled and built.** It does not need re-deriving or
re-proposing. Long design messages restating it are context, not new work; find
the one or two genuine deltas and implement only those.

### Numbers that supersede the README

| quantity | README | actual (2026-08-27) |
|---|---:|---:|
| Arm B events (mini3) | 183 | **141** |
| Arm B events with TSD | 43 | **33** |
| stage 40 universal elements | 310 | **272** |
| stage 40 `BOTH_ARMS` | 21 | **19** |
| stage 40 `ARM_A_MULTICOPY` | 202 | **198** |
| stage 40 `ARM_B_EMPTY_FILLED` | 87 | **55** |

Current stage 40 output is `stage40_v2/`; `stage40_test/` is the stale run the
README quotes. Arm A is unchanged and still accurate in the README: 1,648
families over 200 genomes, ~1.2 s/genome, IS110 recovered at 93% STRONG with
0/13 events carrying a TSD. E00000 remains the IS1-like anchor — 768 bp, 21
copies in one genome, 124/200 genomes, 9 bp TSD.

### Code changes on 2026-08-27

- **Stage 30 no longer deletes oversized bubbles.** They are parked in
  `armB_large_events.tsv` + `armB_large_inserts.fna` with ids `<cid>.L%06d`, and
  `armB_events.tsv` gained a `parent_large_event_id` column. On mini3 this
  recovered **42 bubbles** (median 15 kb, range 5.5–66 kb) that were previously
  discarded at the size gate — including the bin holding the known cases of an
  IS110 nested inside 33 kb and 36 kb cargo.
- **`architecture` is categorical again** (`asymmetric_candidate` /
  `clean_or_graph_slack`). It had been written as `asymmetric_candidate_+1066`
  with the span baked into the label, which made every asymmetric event its own
  singleton class and broke grouping in the stage 61 architecture test. The
  magnitude lives in `junction_span`, where it always was. mini3 is now 64
  asymmetric / 77 clean.
- **Stage 60 rewritten** from a single recall number into a 2×2 contingency
  table. See below.
- **`slurm/isescan_array.sh` added.**

### Stage 60 is a contingency table, not a recall number

ISEScan is not truth — it finds IS by transposase pHMM plus terminal repeats and
therefore carries exactly the family bias this pipeline exists to avoid.

```
                 ISEScan +              ISEScan -
  Arm A +   known IS recovered    ANNOTATION-FREE CANDIDATES  <- the point
  Arm A -   sensitivity failure   background (not enumerable, reported "-")
```

Three design decisions, each load-bearing:

- Genomes lacking either side are dropped. Counting an Arm A family from a
  genome ISEScan never saw would inflate the ArmA+/ISEScan− cell.
- ArmA−/ISEScan+ splits into multi-copy (a real sensitivity failure) and
  singleton (a declared blind spot — Arm A needs ≥3 copies at distinct loci and
  never claimed singletons).
- It does **not** grade its own novel cell. Deciding whether an
  ArmA+/ISEScan− family is real needs the stage 50 screens (barrnap rRNA, ORF
  fraction, transposase HMM), and doing that inside the benchmark would put
  annotation back in front of discovery. Stage 60 isolates and describes the
  pool into `_novel_pool.tsv` and stops.

ArmA+/ISEScan− is **expected to be enriched, not wrong**: IS110/IS1111 has
neither TIR nor TSD — the two features ISEScan keys on — so an annotation-free
method should win there. Reading that cell as false positives inverts the
result.

### Operational lesson: long jobs must go through `sbatch`

Backgrounded shell commands are tied to the session's process lifetime.
`setsid nohup ... &` does **not** survive, and neither do already-running
children. On 2026-08-27 an ISEScan run over 20 genomes was launched as a
background shell, ran ~75 minutes, and was killed at session teardown with
**zero** output — all 20 genomes had started, all 20 were lost.

Anything longer than a few minutes gets a `slurm/*_array.sh` in the repo
convention (`--account=pc_rubinlab --partition=lr6 --qos=lr_normal`,
`FNA_LIST`/`OUT` via `--export`, strided `awk 'NR%n==i'`, and a
skip-if-output-exists guard so a requeue resumes).

ISEScan specifically drives `phmmer --max`, which disables the heuristic filters:
tens of minutes per genome, and it barely scales with `--cpu` — measured at
0.3–0.75 of a core per process on an otherwise 93%-idle 40-core node.
**Parallelise across genomes, not threads.** It also drops scratch
`<dir>_<name>.list` files into the working directory, so each array task needs
its own CWD or concurrent tasks collide.

### Open items, in priority order

1. **The 70 read-validated junctions have never been located.** Stage 61 has
   only ever run against a synthetic gold in `bench61_synth/`. This blocks the
   FNA-first decision and is the one input only the user can supply. The script
   itself is complete: T1=24 / T2=5 / T3=41 tiering, `max(|Δ_L|,|Δ_R|) ≤ k` as
   the headline metric with one-sided agreement reported but explicitly not a
   pass, a separate asymmetric-architecture test, and the four acceptance
   thresholds pre-registered in its docstring.
2. **Stage 32 — the nested sub-5 kb search inside the retained large bubbles —
   is designed but not written.** This is what makes keeping the 42 bubbles pay
   off.
3. **Arm B `cl0000` (9 genomes) returns 0 bubbles** while mini3 (3 genomes)
   returns 248, and its `graph.gfa` is missing from the output dir. Undiagnosed.
4. **Stage 50 has never been run.** No annotation output exists anywhere, and it
   is what grades the stage 60 novel pool.

### In flight at time of writing

SLURM job **25306752**, a 20-task ISEScan array into
`/global/scratch/users/kh36969/fna_ins_discovery/isescan_run`, normalising
per-genome tables to `isescan_run/tsv/<sample>.tsv` for stage 60's
`--isescan-glob`. Check with `squeue -u kh36969`.

### Paths

- Code: `/global/home/users/kh36969/fna_based_mgefinder_project` — **not a git
  repo**, so there is no history to recover from.
- Workdir / outputs: `/global/scratch/users/kh36969/fna_ins_discovery`
- Test genomes: `/global/scratch/users/kh36969/cross_ref_is_runs/results/escherichia_coli/_ncbi_symlinks/*.fna`
  (200 staged *E. coli*)

The PATH that makes stages 10/20/30/30b/40 run with no conda solve:

```bash
export PATH=$PWD/env/bin:$HOME/.conda/envs/claude-env/bin:$PATH
```

`env/bin/` holds static binaries (minigraph 0.21-r606, gfatools 0.5-r296,
skani 0.3.2, pangraph 1.4.0); `claude-env` supplies minimap2, mmseqs, samtools,
nucmer, prodigal, hmmsearch. The `fnains` core env from `setup_env.sh` was never
built and is not needed. Built and working: `fnains_sv` (svim-asm),
`fnains_annot` (isescan.py, barrnap, hmmsearch, seqkit), `fnains_util` (nucdiff,
seqkit).
