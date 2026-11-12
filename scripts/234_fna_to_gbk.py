#!/usr/bin/env python3
"""Convert an assembly FNA to GenBank so GenBank-only tools (ISCompare) can read it.

WHY CDS FEATURES ARE WRITTEN.  ISCompare walks `record.features` of type CDS to
say which gene an insertion lands in.  A featureless GenBank stub parses fine
but makes every gene-context field empty, which would let us conclude the tool
"does not work on our data" when in fact we handed it an empty annotation.  So
genes are predicted with pyrodigal in the same meta mode and translation table
the CDS/Pfam pipeline uses, and written as real CDS features.  --no-genes keeps
the stub behaviour for a deliberate comparison.

The nucleotide sequence is copied verbatim; the conversion adds annotation, it
never edits a base.  Record ids are kept exactly as in the FNA so coordinates
replay against the original file.
"""
import argparse, sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fna", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--organism", default="Bacteria")
    ap.add_argument("--table", type=int, default=11)
    ap.add_argument("--no-genes", action="store_true")
    a = ap.parse_args()

    gf = None
    if not a.no_genes:
        import pyrodigal
        bins = pyrodigal.MetagenomicBins(
            [b for b in pyrodigal.METAGENOMIC_BINS
             if b.training_info.translation_table == a.table])
        gf = pyrodigal.GeneFinder(meta=True, metagenomic_bins=bins,
                                  closed=False, mask=False)

    out, ncds, nrec, nbp = [], 0, 0, 0
    for rec in SeqIO.parse(a.fna, "fasta"):
        s = str(rec.seq).upper()
        r = SeqRecord(Seq(s), id=rec.id, name=rec.id[:16],
                      description=rec.description or rec.id)
        r.annotations["molecule_type"] = "DNA"
        r.annotations["organism"] = a.organism
        r.annotations["source"] = a.organism
        r.annotations["topology"] = "linear"      # we do not know; recorded as unknown upstream
        if gf is not None and len(s) >= 20:
            for i, g in enumerate(gf.find_genes(s.encode()), 1):
                f = SeqFeature(
                    FeatureLocation(g.begin - 1, g.end, strand=g.strand),
                    type="CDS",
                    qualifiers={"locus_tag": ["%s_%05d" % (rec.id.replace("|", "_"), i)],
                                "product": ["hypothetical protein"],
                                "transl_table": [str(a.table)],
                                "translation": [g.translate(include_stop=False)]})
                r.features.append(f); ncds += 1
        out.append(r); nrec += 1; nbp += len(s)
    SeqIO.write(out, a.out, "genbank")
    sys.stderr.write("%s: %d records, %d bp, %d CDS -> %s\n"
                     % (a.fna.split("/")[-1], nrec, nbp, ncds, a.out))

main()
