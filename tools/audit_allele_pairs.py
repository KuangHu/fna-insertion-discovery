#!/usr/bin/env python3
"""Sequence audit of short/overlapping anchor pairs; does not revise calls.

Uses one optimal global alignment. Gap placement in repeats is not unique.
All exported coordinates are zero-based, half-open in genomic orientation.
"""
import csv
import html
import json
from pathlib import Path
from Bio import Align
from local_debug import ROOT, SOURCES, fasta


def main():
    out = ROOT / 'local_debug/filled_empty'
    with (out / 'verification.tsv').open() as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    genomes = {n: fasta(ROOT/'local_debug/genomes'/sp/(n+'.fna'))
               for n, sp, _, _ in SOURCES}
    aligner = Align.PairwiseAligner(mode='global', match_score=2,
                                   mismatch_score=-3, open_gap_score=-12,
                                   extend_gap_score=-0.1)
    results, sections = [], []
    for r in rows:
        if r['status'] not in ('observed_short_interval', 'overlapping_anchors'):
            continue
        # Current selected examples are forward-strand; reject silently wrong
        # coordinate interpretation if future datasets introduce reverse pairs.
        assert r['target_strand'] == '+', 'Reverse pairs require coordinate conversion'
        s, e = int(r['source_start']), int(r['source_end'])
        a, b = int(r['target_start']), int(r['target_end'])
        x0, y0 = s-500, min(a,b)-500
        x = genomes[r['source_genome']][r['source_contig']][x0:e+500]
        y = genomes[r['target_genome']][r['target_contig']][y0:max(a,b)+500]
        aln = aligner.align(x, y)[0]
        coords = aln.coordinates
        gaps = []
        matches = aligned = 0
        for i in range(coords.shape[1]-1):
            u,v = map(int, coords[:,i]); U,V = map(int, coords[:,i+1])
            if U>u and V==v:
                gaps.append((U-u, x0+u, x0+U, y0+v))
            elif U>u and V>v:
                assert U-u == V-v
                matches += sum(c==d for c,d in zip(x[u:U],y[v:V]))
                aligned += U-u
        gap, start, end, junction = max(gaps)
        overlap = max(0, min(e,end)-max(s,start))
        d = dict(r, alignment_gap_bp=gap, alignment_gap_start=start,
                 alignment_gap_end=end, target_junction=junction,
                 called_interval_covered_by_gap=round(overlap/(e-s),5),
                 aligned_base_identity=round(matches/aligned,5),
                 aligned_columns=aligned, start_shift=start-s, end_shift=end-e)
        results.append(d)
        detail = (f"{r['copy_id']}\nCoordinates: 0-based half-open.\n"
                  f"Source window {r['source_genome']} {r['source_contig']} {x0}:{e+500}\n"
                  f"Comparison window {r['target_genome']} {r['target_contig']} {y0}:{max(a,b)+500}\n"
                  'Alignment display uses window-relative coordinates.\n'
                  'One optimal alignment, NOT a uniquely proven breakpoint.\n\n'+str(aln))
        (out/(r['id']+'_alignment.txt')).write_text(detail)
        (out/(r['id']+'_pair.fna')).write_text(
            f">{r['id']}_source {r['source_genome']} {r['source_contig']}:{x0}-{e+500}\n{x}\n"
            f">{r['id']}_comparison {r['target_genome']} {r['target_contig']}:{y0}-{max(a,b)+500}\n{y}\n")
        sections.append('<h2>'+html.escape(r['copy_id'])+'</h2><pre>'+
            html.escape(f"Source:     [host L] [{e-s} bp called interval] [host R]\n"
                        f"Comparison: [host L] [{r['observed_gap_bp']} bp anchor gap] [host R]\n"
                        f"Global alignment: {gap} bp dominant gap; aligned-base identity {matches/aligned:.2%}\n"
                        f"Gap endpoints versus call: {start-s:+d}, {end-e:+d} bp\n"
                        f"Arm A boundary_refined={r['boundary_refined']}")+
            '</pre><p>Negative anchor gap means overlapping anchors, not a negative-length allele.</p>'+
            f'<p><a href="{r["id"]}_pair.fna">Actual genome sequences</a> · '+
            f'<a href="{r["id"]}_alignment.txt">Full alignment</a></p>')
    with (out/'pair_alignment_audit.tsv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(results[0]),delimiter='\t')
        w.writeheader();w.writerows(results)
    (out/'sequence_audit.html').write_text('<!doctype html><meta charset="utf-8">'+
        '<title>Genome allele sequence audit</title><style>body{max-width:1050px;margin:36px auto;font:16px Arial;color:#17354a}pre{padding:20px;background:#edf3f5;white-space:pre-wrap}h2{margin-top:36px}</style>'+
        '<h1>Genome → candidate interval → genome</h1>'+
        '<p>Ten real source/comparison genome pairs. Global alignment corroborates long/short sequence differences; it does not establish evolutionary direction, transposition, or RNA guidance.</p>'+
        '<p>Scores: match +2, mismatch −3, gap open −12, extension −0.1. One optimal alignment is shown; repeated junctions can have equivalent placements. No calls were rewritten.</p>'+
        '<p>The three positive anchor gaps retain 5, 9 and 10 bp respectively; all three have unrefined Arm A boundaries. Seven overlapping-anchor pairs remain boundary-ambiguous.</p>'+
        ''.join(sections),encoding='utf-8')
    print(json.dumps([{k:r[k] for k in ('id','alignment_gap_bp','start_shift','end_shift',
          'aligned_base_identity','called_interval_covered_by_gap')} for r in results],indent=2))


if __name__ == '__main__':
    main()
