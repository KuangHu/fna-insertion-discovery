#!/usr/bin/env python3
"""Audit every Arm B event against actual carrier genomes, including large events."""
import csv
import hashlib
import html
import importlib.util
import json
from pathlib import Path
import sys
from Bio import Align, SeqIO
from classify_allele_evidence import compare, rc
from local_debug import ROOT

sys.path.insert(0,str(ROOT/'scripts'))
spec=importlib.util.spec_from_file_location('armb',ROOT/'scripts/30_armB_graph.py')
armb=importlib.util.module_from_spec(spec);spec.loader.exec_module(armb)


def read(path):
    with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))


def main():
    root=ROOT/'local_debug/armB'
    genomes={n:{r.id:str(r.seq) for r in SeqIO.parse(ROOT/'local_debug/genomes/spneumoniae'/(n+'.fna'),'fasta')}
             for n in ('spneumoniae_R6','spneumoniae_TIGR4')}
    allrows=[];summary={}
    local=Align.PairwiseAligner(mode='local',match_score=2,mismatch_score=-3,open_gap_score=-12,extend_gap_score=-.1)
    for panel in ('r6_first','tigr4_first'):
        folder=root/panel;out=folder/'audit';out.mkdir(exist_ok=True)
        calls={n:armb.parse_call(folder/'calls'/(n+'.bed')) for n in genomes}
        for category,table,fa in [('primary','armB_events.tsv','armB_inserts.fna'),
                                  ('large','armB_large_events.tsv','armB_large_inserts.fna')]:
            inserts={r.id:str(r.seq) for r in SeqIO.parse(folder/fa,'fasta')}
            for event in read(folder/table):
                eid=event.get('event_id',event.get('large_event_id'))
                key=(event['chrom'],int(event['bubble_start']),int(event['bubble_end']))
                fn=event['filled_genomes'].split(',')[0];en=event['empty_genomes'].split(',')[0]
                fl,fc,fs,fe,fstrand,_=calls[fn][key]
                el,ec,es,ee,estrand,_=calls[en][key]
                fs,fe=sorted((fs,fe));es,ee=sorted((es,ee))
                fseq=genomes[fn][fc];eseq=genomes[en][ec]
                assert 0<=fs<=fe<=len(fseq) and 0<=es<=ee<=len(eseq)
                x=fseq[fs:fe];y=eseq[es:ee]
                graphinsert=inserts[eid]
                assert hashlib.md5(graphinsert.upper().encode()).hexdigest()==event['insert_md5']
                oriented=x if fstrand=='+' else rc(x)
                matches=columns=covered=0
                if graphinsert in oriented:
                    matches=columns=covered=len(graphinsert)
                else:
                    aln=local.align(oriented,graphinsert)[0]
                    for j in range(aln.coordinates.shape[1]-1):
                        a,b=map(int,aln.coordinates[:,j]);A,B=map(int,aln.coordinates[:,j+1])
                        if A>a and B>b:
                            columns+=A-a;covered+=B-b
                            matches+=sum(u==v for u,v in zip(oriented[a:A],graphinsert[b:B]))
                row=dict(event_id=eid,panel=panel,category=category,filled_genome=fn,filled_contig=fc,
                         filled_start=fs,filled_end=fe,short_genome=en,short_contig=ec,short_start=es,short_end=ee,
                         graph_length_difference=event['insert_size'],genome_interval_difference=len(x)-len(y),
                         filled_interval_bp=len(x),short_interval_bp=len(y),exported_fragment_bp=len(graphinsert),
                         fragment_identity=round(matches/columns,6),fragment_coverage=round(covered/len(graphinsert),6),
                         sequence_support='supported' if matches/columns>=.95 and covered/len(graphinsert)>=.95 else 'needs_review',
                         long_short_evidence='not_aligned_large',gap_bp='.',gap_call_coverage='.',aligned_identity='.')
                x0,x1=max(0,fs-500),min(len(fseq),fe+500)
                y0,y1=max(0,es-500),min(len(eseq),ee+500)
                X,Y=fseq[x0:x1],eseq[y0:y1]
                reverse=fstrand!=estrand
                if reverse:Y=rc(Y)
                (out/(eid+'_pair.fna')).write_text(
                    f">long {fn} {fc}:{x0}-{x1} strand=+\n{X}\n>short {en} {ec}:{y0}-{y1} strand={'-' if reverse else '+'}\n{Y}\n")
                if category=='primary':
                    metrics,alignment=compare(X,Y,fs-x0,fe-x0)
                    row.update(long_short_evidence=metrics['evidence_status'],gap_bp=metrics['dominant_gap_bp'],
                               gap_call_coverage=metrics['gap_call_coverage'],aligned_identity=metrics['aligned_base_identity'])
                    (out/(eid+'_alignment.txt')).write_text('Window-relative coordinates. One optimal alignment, not a proven mechanism.\n'+alignment)
                allrows.append(row)
        from collections import Counter
        panelrows=[r for r in allrows if r['panel']==panel]
        summary[panel]=dict(bubbles=len(read(folder/'bubbles.tsv')),categories=dict(Counter(r['category'] for r in panelrows)),
                            sequence_support=dict(Counter(r['sequence_support'] for r in panelrows)),
                            primary_alignment=dict(Counter(r['long_short_evidence'] for r in panelrows if r['category']=='primary')))
    with (root/'audit.tsv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(allrows[0]),delimiter='\t');w.writeheader();w.writerows(allrows)
    # Compare loci in R6 coordinates, regardless of whether R6 is the long carrier.
    def interval(r):
        return (r['filled_start'],r['filled_end']) if r['filled_genome']=='spneumoniae_R6' else (r['short_start'],r['short_end'])
    first=[r for r in allrows if r['panel']=='r6_first' and r['category']=='primary']
    second=[r for r in allrows if r['panel']=='tigr4_first' and r['category']=='primary']
    edges=[]
    for a in first:
        s,e=interval(a)
        for b in second:
            u,v=interval(b)
            if a['filled_genome']==b['filled_genome'] and abs(s-u)<=100 and abs(e-v)<=100:
                edges.append((a['event_id'],b['event_id']))
    summary['order_comparison']=dict(rule='Same long carrier, both R6 genomic endpoints within 100 bp; proximity comparison only',
        paired_edges=len(edges),r6_first_matched=len(set(a for a,b in edges)),tigr4_first_matched=len(set(b for a,b in edges)))
    ea=read(root/'r6_first/armB_events.tsv');eb=read(root/'tigr4_first/armB_events.tsv')
    summary['order_comparison']['identical_primary_fragment_md5_multiset']=(Counter(r['insert_md5'] for r in ea)==Counter(r['insert_md5'] for r in eb))
    summary['order_comparison']['same_order_fragment_md5_matches']=sum(a['insert_md5']==b['insert_md5'] for a,b in zip(ea,eb))
    for panel in ('r6_first','tigr4_first'):
        rr=[r for r in allrows if r['panel']==panel and r['category']=='primary']
        delta=[abs(int(r['graph_length_difference'])-r['genome_interval_difference']) for r in rr]
        summary[panel]['graph_vs_genome_length_disagreements']=sum(d!=0 for d in delta)
        summary[panel]['maximum_length_difference_bp']=max(delta)
    (root/'order_pairs.tsv').write_text('r6_first_event\ttigr4_first_event\n'+''.join(a+'\t'+b+'\n' for a,b in edges))
    (root/'audit_summary.json').write_text(json.dumps(summary,indent=2))
    cols=list(allrows[0])
    body=''.join('<tr>'+''.join('<td>'+html.escape(str(r[k]))+'</td>' for k in cols)+
                 f'<td><a href="{r["panel"]}/audit/{r["event_id"]}_pair.fna">Genome pair</a></td></tr>' for r in allrows)
    (root/'audit.html').write_text('<!doctype html><meta charset="utf-8"><title>Arm B genome audit</title><style>body{font:15px Arial;margin:30px}td,th{padding:8px;border-bottom:1px solid #ddd}table{border-collapse:collapse}pre{white-space:pre-wrap}</style><h1>Arm B: original genome audit</h1><pre>'+html.escape(json.dumps(summary,indent=2))+'</pre><p>All primary events aligned with 500 bp flanks. Large events have coordinate, sequence-export and local sequence checks only. Bubble intervals may include shared sequence; a failed dominant-gap threshold is not absence. Fragment identity excludes gaps; coverage counts aligned fragment bases. This is assembly-level support, not transposition or RNA-guidance validation.</p><table><tr>'+''.join('<th>'+k+'</th>' for k in cols)+'<th>Evidence</th></tr>'+body+'</table>',encoding='utf-8')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
