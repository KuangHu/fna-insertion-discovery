#!/usr/bin/env python3
"""Apply identical allele evidence rules to every mapped candidate copy.

Consumes verify_filled_empty.py output. Does not change biological call boundaries.
"""
import collections
import argparse
import csv
import html
import json
from pathlib import Path
from Bio import Align, SeqIO
from local_debug import ROOT, SOURCES


def rc(seq):
    return seq.translate(str.maketrans('ACGTRYMKSWBDHVNacgtrymkswbdhvn',
                                     'TGCAYRKMSWVHDBNtgcayrkmswvhdbn'))[::-1]


def classify(coverage, identity, left_columns, right_columns):
    # Operational evidence thresholds, not calibrated biological probabilities.
    return ('supported_long_short' if coverage >= .90 and identity >= .95
            and min(left_columns, right_columns) >= 200
            else 'no_dominant_long_short_support')


def compare(x, y, called_start, called_end):
    if not (0 <= called_start < called_end <= len(x)) or not y:
        raise ValueError('Invalid called interval or empty comparison window')
    x,y=x.upper(),y.upper()
    aligner = Align.PairwiseAligner(mode='global', match_score=2,
                                   mismatch_score=-3, open_gap_score=-12,
                                   extend_gap_score=-.1)
    aln = aligner.align(x, y)[0]
    gaps = []
    matches = columns = left = right = 0
    for i in range(aln.coordinates.shape[1]-1):
        u,v = map(int, aln.coordinates[:,i])
        U,V = map(int, aln.coordinates[:,i+1])
        if U>u and V==v:
            overlap = max(0, min(U,called_end)-max(u,called_start))
            gaps.append((overlap,U-u,u,U))
        elif U>u and V>v:
            assert U-u == V-v
            matches += sum(a==b for a,b in zip(x[u:U],y[v:V]))
            columns += U-u
            left += max(0,min(U,called_start)-u)
            right += max(0,U-max(u,called_end))
    overlap,gap,start,end = max(gaps,default=(0,0,called_start,called_start))
    coverage = overlap/(called_end-called_start)
    identity = matches/columns if columns else 0
    return dict(evidence_status=classify(coverage,identity,left,right),
                dominant_gap_bp=gap, gap_call_coverage=round(coverage,6),
                aligned_base_identity=round(identity,6),
                left_aligned_bp=left,right_aligned_bp=right,
                gap_start_shift=start-called_start if gap else '.',
                gap_end_shift=end-called_end if gap else '.'), str(aln)


def validate(rows, genomes):
    if not rows:
        raise ValueError('Verification table has no rows')
    ids=set()
    for r in rows:
        ident=r['id']
        if not ident or any(c not in 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-' for c in ident):
            raise ValueError('id must be a filename-safe unique row identifier')
        if ident in ids:
            raise ValueError('Duplicate row id: '+ident)
        ids.add(ident)
        source=genomes[r['source_genome']][r['source_contig']]
        s,e=int(r['source_start']),int(r['source_end'])
        if not 0<=s<e<=len(source) or int(r['source_insert_bp'])!=e-s:
            raise ValueError('Invalid source coordinates/length: '+ident)
        if r['boundary_refined'] not in ('0','1'):
            raise ValueError('Invalid boundary_refined: '+ident)
        count=int(r['qualifying_pairs'])
        if count<0 or (r['status']=='no_comparator' and count):
            raise ValueError('Invalid pair count: '+ident)
        if count==1:
            target=genomes[r['target_genome']][r['target_contig']]
            a,b=int(r['target_start']),int(r['target_end'])
            if not (0<=a<=len(target) and 0<=b<=len(target)):
                raise ValueError('Invalid comparison coordinates: '+ident)
            if r['target_strand'] not in ('+','-'):
                raise ValueError('Invalid comparison strand: '+ident)
            if int(r['observed_gap_bp'])!=b-a:
                raise ValueError('Inconsistent anchor gap: '+ident)


def run(rows, genomes, out):
    """Write the same evidence products for all rows after validating inputs."""
    validate(rows,genomes)
    out=Path(out)
    out.mkdir(parents=True,exist_ok=True)
    results=[]
    with (out/'uniform_middle_sequences.fna').open('w') as middle:
        for r in rows:
            s,e=int(r['source_start']),int(r['source_end'])
            source=genomes[r['source_genome']][r['source_contig']]
            middle.write(f">{r['id']}_called_interval\n{source[s:e]}\n")
            d=dict(r,comparison_interval_bp='.',comparison_interval_sequence='.',
                   evidence_status='unresolved',boundary_status='not_assessable',
                   dominant_gap_bp='.',gap_call_coverage='.',aligned_base_identity='.',
                   left_aligned_bp='.',right_aligned_bp='.',gap_start_shift='.',gap_end_shift='.')
            if r['status']=='no_comparator':
                d['evidence_status']='no_comparator'
            elif int(r['qualifying_pairs'])==1:
                a,b=int(r['target_start']),int(r['target_end'])
                target=genomes[r['target_genome']][r['target_contig']]
                reverse=r['target_strand']=='-'
                if b>=a:
                    seq=target[a:b]
                    if reverse:seq=rc(seq)
                    d.update(comparison_interval_bp=b-a,
                             comparison_interval_sequence=seq or '(zero-length)')
                    if seq:middle.write(f">{r['id']}_comparison_interval\n{seq}\n")
                x0,x1=max(0,s-500),min(len(source),e+500)
                y0,y1=max(0,min(a,b)-500),min(len(target),max(a,b)+500)
                x,y=source[x0:x1],target[y0:y1]
                if reverse:y=rc(y)
                metrics,alignment=compare(x,y,s-x0,e-x0)
                d.update(metrics)
                d['boundary_status']=('anchor_overlap' if b<a else
                    'nominal_boundary' if r['boundary_refined']=='0' else
                    'refined_boundary_not_mechanistically_verified')
                (out/(r['id']+'_uniform_alignment.txt')).write_text(
                    f"Source {r['source_genome']} {r['source_contig']} [{x0},{x1}) +\n"
                    f"Comparison {r['target_genome']} {r['target_contig']} [{y0},{y1}) {r['target_strand']}\n"
                    'Display coordinates are relative to oriented windows; one optimal alignment.\n'+alignment)
            results.append(d)
    with (out/'uniform_evidence.tsv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(results[0]),delimiter='\t')
        w.writeheader();w.writerows(results)
    summary=dict(total=len(results),counts=dict(collections.Counter(r['evidence_status'] for r in results)),
        policy='Same rules for every copy. >=90% call coverage by one source-only gap, >=95% aligned-base identity, >=200 aligned bp in each source flank. Operational thresholds, not validated probabilities.',
        limits='Short sequences preserved; overlaps are undefined middle intervals, never negative alleles. Unresolved is not absent. No mechanistic or evolutionary polarity claims.')
    (out/'uniform_summary.json').write_text(json.dumps(summary,indent=2))
    cols=['copy_id','evidence_status','boundary_status','source_insert_bp','comparison_interval_bp',
          'dominant_gap_bp','gap_start_shift','gap_end_shift']
    body=''.join('<tr>'+''.join('<td>'+html.escape(str(r[k]))+'</td>' for k in cols)+'</tr>' for r in results)
    (out/'uniform_report.html').write_text('<!doctype html><meta charset="utf-8"><title>Uniform allele evidence</title>'+
        '<style>body{font:16px Arial;margin:32px;color:#17354a}table{border-collapse:collapse}td,th{padding:9px;border-bottom:1px solid #ddd;text-align:left}pre{white-space:pre-wrap}</style>'+
        '<h1>Uniform allele evidence</h1><p>[Host L] [called interval] [Host R] versus [homologous L] [observed interval] [homologous R]</p>'+
        '<p>Evidence and boundary status are independent. TTTTT and ATTATATAT are preserved as observed sequence, without motif or TSD labels.</p><pre>'+
        html.escape(json.dumps(summary,indent=2))+'</pre><p><a href="uniform_evidence.tsv">Full table</a> · <a href="uniform_middle_sequences.fna">Middle sequences</a></p><table><tr>'+
        ''.join('<th>'+k+'</th>' for k in cols)+'</tr>'+body+'</table>',encoding='utf-8')
    print(json.dumps(summary,indent=2))
    return summary


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verification',type=Path,help='Anchor-pair verification TSV')
    parser.add_argument('--manifest',type=Path,help='TSV: genome_id, fasta_path; paths relative to manifest')
    parser.add_argument('--outdir',type=Path,help='Output directory')
    args=parser.parse_args()
    if any((args.verification,args.manifest,args.outdir)) and not all((args.verification,args.manifest,args.outdir)):
        parser.error('Provide --verification, --manifest and --outdir together')
    out=args.outdir or ROOT/'local_debug/filled_empty'
    with (args.verification or out/'verification.tsv').open() as f:
        rows=list(csv.DictReader(f,delimiter='\t'))
    paths={}
    if args.manifest:
        with args.manifest.open() as f:
            for r in csv.DictReader(f,delimiter='\t'):
                name=r['genome_id']
                if name in paths:raise ValueError('Duplicate genome_id: '+name)
                paths[name]=args.manifest.parent/Path(r['fasta_path'])
    else:
        paths={n:ROOT/'local_debug/genomes'/sp/(n+'.fna') for n,sp,_,_ in SOURCES}
    genomes={n:{k:str(v.seq) for k,v in SeqIO.to_dict(SeqIO.parse(p,'fasta')).items()}
             for n,p in paths.items()}
    run(rows,genomes,out)


if __name__=='__main__':main()
