#!/usr/bin/env python3
"""Concatenate the genome corpus under serial ids, recording contig lengths.

The length is not cosmetic: the context extractor needs it to know whether a
2000bp flank fits before the contig ends, and a hit near an end must be
reported as truncated rather than silently short.
"""
import os
import sys

lst, out_fa, out_map = sys.argv[1], sys.argv[2], sys.argv[3]
n = 0
with open(out_fa, "w") as fa, open(out_map, "w") as mp:
    for gi, path in enumerate(l.strip() for l in open(lst) if l.strip()):
        genome = os.path.basename(path)
        for ext in (".fna", ".fa", ".fasta"):
            if genome.endswith(ext):
                genome = genome[:-len(ext)]
                break
        contig, parts = None, []

        def flush():
            global n
            if contig is None:
                return
            seq = "".join(parts)
            n += 1
            sid = "S%08d" % n
            fa.write(">%s\n" % sid)
            for i in range(0, len(seq), 80):
                fa.write(seq[i:i + 80] + "\n")
            mp.write("%s\t%s\t%s\t%d\n" % (sid, genome, contig, len(seq)))

        try:
            fh = open(path)
        except OSError as e:
            print("WARN: cannot open %s (%s)" % (path, e), file=sys.stderr)
            continue
        with fh:
            for line in fh:
                if line.startswith(">"):
                    flush()
                    contig = line[1:].split()[0]
                    parts = []
                else:
                    parts.append(line.strip())
            flush()
            contig = None
        if (gi + 1) % 2000 == 0:
            print("  %d genomes, %d contigs" % (gi + 1, n), file=sys.stderr)
print("corpus: %d contigs" % n, file=sys.stderr)
