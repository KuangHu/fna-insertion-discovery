# Lustre scratch I/O errors on cf1 nodes, 2026-11-10 to 2026-11-11

Brief for LRC support. Account `pc_rubinlab`, partition `cf1`, QOS `cf_normal`, user `kh36969`.
All paths are under `/global/scratch/users/kh36969/fna_ins_discovery/isref_all/`, written as `$W/` below.
Times are PST (UTC−8). Error lines are quoted verbatim from the job logs.

## Summary

- Two jobs on nine cf1 nodes (11-10/11; more on 11-11 evening, see the update section) got `Input/output error` (EIO, errno 5) or `Cannot send after transport endpoint shutdown` (ESHUTDOWN, errno 108). This happened on both reads and writes to `/global/scratch`.
- The errors are intermittent. They hit several unrelated nodes in the same evening (roughly 20:30–02:00 on 11-10/11), and one node again on 11-11 afternoon.
- On n0052 (job 27786542), every read of six existing files failed partway through.
- One task did not end after a time-limit cancel: slurmstepd reported `JOB NOT ENDING WITH SIGNALS`. A file in that task's output path has an mtime 2 h 14 min after the job ended.
- We have also seen one possible silent corruption: a line in a file written on scratch had two records merged into one, and no error was reported for that write (see job 27725308, task 8, below).

## Questions for LRC

1. Were there Lustre client evictions or OST/MDT problems on these nodes in these windows?
   - n0037, n0039, n0043, n0044, n0046, n0047 and n0051, 11-10 20:00 – 11-11 02:10;
   - n0040 around 11-11 02:00–04:20;
   - n0052, 11-11 17:08–17:35.
2. Are OST 9–11 known to be slow (3–5× slower per file in our probe)? Is n0052 still unhealthy? Should we keep excluding it, and are there other nodes we should exclude?
3. Task 4 of job 27725308 (n0040): processes were left behind after `JOB NOT ENDING WITH SIGNALS`. Could they have kept writing to scratch after the job ended?

## Job 27786542 (single job, node n0052.cf1, 11-11 17:08:02 – 17:34:35, FAILED, exit 123)

- What it did: a read-only `sha256sum` of six existing files (24 GB in total), run in parallel with `xargs -P 6`.
- All six reads failed partway through. The log (stderr lines interleaved by the parallel processes; unmodified) reads:

```
sha256sum: hits/shard_50.rawsha256sum: hits/shard_55.tsv.gzsha256sum: copies/shard_50.tsvsha256sum: hits/shard_52.rawsha256sum: hits/shard_55.tsv: Input/output error: Input/output error: Cannot send after transport endpoint shutdown: Cannot send after transport endpoint shutdown
: Input/output error



sha256sum: hits/shard_50.tsv: Input/output error
```

- The files are unchanged: their size, mtime and mode were the same before and after the job.

## Job 27725308 (array of 16 tasks, 11-10 12:00 – 11-11 02:00)

Only tasks with I/O-related lines are listed. Each "approx. time" is the mtime of the last file written before the error, so it is only approximate.

| Task (SLURM job id) | Node | Approx. time | Verbatim error |
|---|---|---|---|
| 1 (27725314) | n0037.cf1 | before 20:23 (task end) | `/var/spool/slurmd/job27725314/slurm_script: line 53: $W/shardlog/FAILED_33: Cannot send after transport endpoint shutdown` |
| 3 (27725316) | n0039.cf1 | 22:49 | `awk: cmd. line:1: (FILENAME=$W/hits/shard_35.raw FNR=21834778) fatal: print to "standard output" failed (Input/output error)` |
| 3 (27725316) | n0039.cf1 | 23:20 | `OSError: [Errno 5] Input/output error` (Python, reading `hits/shard_35.tsv`) |
| 3 (27725316) | n0039.cf1 | after 23:20, before 02:00 | 18 lines: `/var/spool/slurmd/job27725316/slurm_script: line 53: $W/shardlog/FAILED_51: Cannot send after transport endpoint shutdown` (16 times) and `… FAILED_51: Input/output error` (2 times) |
| 7 (27725320) | n0043.cf1 | 01:57:54 | `gzip: $W/hits/shard_55.tsv: Input/output error` |
| 10 (27725323) | n0046.cf1 | 21:27 | `OSError: [Errno 5] Input/output error` (Python, reading `hits/shard_42.tsv`) |
| 10 (27725323) | n0046.cf1 | 23:47 | `awk: cmd. line:1: (FILENAME=$W/hits/shard_58.raw FNR=1453640) fatal: print to "standard output" failed (Input/output error)` |
| 11 (27725324) | n0047.cf1 | before 21:44 | `awk: cmd. line:1: (FILENAME=$W/hits/shard_43.raw FNR=1389289) fatal: print to "standard output" failed (Input/output error)` and `/var/spool/slurmd/job27725324/slurm_script: line 75: $W/shardlog/chain_43.log: Input/output error` |
| 11 (27725324) | n0047.cf1 | about 21:44 | `/var/spool/slurmd/job27725324/slurm_script: line 72: $W/hits/shard_59.tsv: Input/output error` |
| 15 (array task id 27725308) | n0051.cf1 | 22:16 | `BrokenPipeError: [Errno 108] Cannot send after transport endpoint shutdown` (Python, `open(... "w")` of a new file in `copies/`) |

Notes on the table:
- **Task 1 (n0037).** Line 53 of the script appends a line to a log file after a genome FASTA under `/global/scratch` could not be read. The same task reported two shards as `632/633 genomes`, so on this node reads failed for at least two input files as well as the log writes.
- **Task 3 (n0039).** The 18 failed appends on line 53 likewise mean that 18 genome-file reads failed, followed by 18 failed log writes.
- **Task 7 (n0043).** gzip stopped partway: `hits/shard_55.tsv.gz` is 728,236,032 bytes and still has gzip's temporary mode `-rw-------`.
- **Task 11 (n0047).** `hits/shard_59.raw` stopped at exactly 553,648,128 bytes (528 MiB). On n0039, `hits/shard_35.tsv` stopped at exactly 2,147,483,648 bytes (2^31). Both are suspiciously round sizes.

**Task 8 (27725321, n0044.cf1): possible silent corruption.**
- A 5.3 GB file written on scratch (`hits/shard_56.tsv`, written by `awk`, mtime 23:11) contains a line in which the end of one record and the start of the next are missing. Two IDs are fused into `Q_5a87aa5Q_5b08ecd47904615a`.
- No error was reported for that write.
- We cannot tell whether this happened in the filesystem or in our pipeline. We report it because it falls in the same window.

**Task 4 (27725317, n0040.cf1): job did not end after the time-limit cancel.**

```
[2026-11-11T02:00:54.503] error: *** JOB 27725317 ON n0040.cf1 CANCELLED AT 2026-11-11T02:00:54 DUE TO TIME LIMIT ***
[2026-11-11T02:05:55.000] error: *** JOB 27725317 STEPD TERMINATED ON n0040.cf1 AT 2026-11-11T02:05:55 DUE TO JOB NOT ENDING WITH SIGNALS ***
[2026-11-11T02:06:13.533] error: Container 66139 in cgroup plugin has 2 processes, giving up after 319 sec
```

- `$W/hits/shard_52.raw` was being written by this task (blastn output). It is 133,572,473 bytes with mtime **2026-11-11 04:15:09**, 2 h 14 min after the job was cancelled.
- None of our jobs were running at that time, so we do not know what wrote to it. Possibly the two leftover processes were stuck in uninterruptible I/O and finished a pending write once the filesystem recovered.

## Update, 2026-11-11 evening

- **Job 27788009** (n0057.cf1, 18:46:32 – 18:48:35): `CANCELLED by 0`, i.e. by root. The batch script never started: its stdout file `$W/v0/logs/isall_V0_27788009.out` was never created. The previous job on that node (27786705, another user) ended FAILED at 18:33.
- **Job 27788052** (n0059.cf1, 18:50:10 – 19:45:04, FAILED): reads genome FASTA files from `/global/scratch` with 64 workers.
  - In 55 minutes it read only 38.4 GB (`MaxDiskRead`). It used 1 h 20 min of CPU time, about 2% of 64 cores, so the workers were mostly blocked on I/O.
  - It then hit 25 read errors of type EIO, ENOTCONN or ESHUTDOWN in a short burst, and our tool stopped at its limit of 5.
  - The verbatim failure record is `node_io_error_limit: node n0059.cf1 had 25 I/O errors (limit 5)`.
  - sacct shows six other FAILED jobs from another user on n0059.cf1 between 16:56 and 18:44 the same day.

## Read probe, 2026-11-11 20:28–20:31 (job 27788912, n0012.cf1)

- **What we ran.** A read-only probe on the same genome files:
  - read 3,240 randomly chosen files in full, with 1, 8, 32 and 64 parallel readers;
  - ran `stat` on all 40,541 files (179.7 GB in total);
  - wrote 512 MiB with `fsync`, then read it back.
- **Errors: none.**
  - Throughput at 64 readers was 1,333 MB/s. Job 27788052 on n0059.cf1 two hours earlier averaged about 12 MB/s.
  - The `stat` pass took 37.5 s, with p99 0.24 s per file.
- **Per-OST read speed.** Median per file, with the 10th percentile in brackets:

  | OSTs | Speed |
  |---|---|
  | OST 9 / 10 / 11 (the 28%-full targets) | 12 / 12.5 / 8 MB/s (4 / 5 / 3) |
  | Every other HDD OST (0–8, 12–21) | 30–61 MB/s (15–33) |

  OST 9, 10 and 11 are 3–5× slower than the rest; 145–165 files were sampled on each.
- **One stall.** Two reads on different OSTs (10 and 12) each took 50.3–50.5 s, ending together at about 20:29:40. All other reads finished in under 1.3 s.
- **Our reading.** Both reads stalled at the same moment on two different OSTs, which looks like a stall on the client or across the filesystem. A longer stall of the same kind would match the EIO/ENOTCONN/ESHUTDOWN errors above. The problem comes and goes over time; we have not found one node or one OST that is permanently broken.

## Our mitigations (for context)

- The next job excludes n0052 (`--exclude=n0052`).
- Our genome reader now re-reads a file up to 3 times after EIO, ENOTCONN or ESHUTDOWN, with increasing waits. It accepts a file only after two complete reads agree byte for byte, and it stops the whole run after 5 such errors on one node.
- The affected outputs are being regenerated. The partial files are kept read-only and are not used.
