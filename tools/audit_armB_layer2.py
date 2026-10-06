#!/usr/bin/env python3
"""Compare locked Layer 2 runs without conflating padded alleles with inserts."""
import collections
import csv
import hashlib
import html
import json
import sys
from Bio import SeqIO
from local_debug import ROOT
sys.path.insert(0,str(ROOT/'scripts'))
from lib_insert import canonical_insert_key, rc


def read(path):
    with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))


def main():
    root=ROOT/'local_debug/armB'
    genomes={n:{r.id:str(r.seq) for r in SeqIO.parse(ROOT/'local_debug/genomes/spneumoniae'/(n+'.fna'),'fasta')}
             for n in ('spneumoniae_R6','spneumoniae_TIGR4')}
    results={};summary={};sequences_checked=0
    for panel in ('r6_first','tigr4_first'):
        folder=root/panel/'layer2'
        rows=read(folder/'loci.tsv');results[panel]={r['locus_id'].split('.')[-1]:r for r in rows}
        seqs={tuple(r.id.split('|')[:2]):str(r.seq) for r in SeqIO.parse(folder/'alleles.fna','fasta')}
        for a in read(folder/'alleles.tsv'):
            seq=seqs[a['locus_id'],a['allele_id']]
            if seq=='-':seq=''
            assert len(seq)==int(a['length'])
            assert hashlib.md5(seq.encode()).hexdigest()[:12]==a['md5']
            for name in a['genomes'].split(','):
                assert any(seq in genome or rc(seq) in genome for genome in genomes[name].values())
                sequences_checked+=1
        for r in rows:
            long=seqs[r['locus_id'],r['longest_allele']]
            start=int(r['insert_start_in_long_bp']);length=int(r['inserted_len'])
            fragment=long[start:start+length]
            digest=hashlib.md5(fragment.encode()).hexdigest()[:12] if fragment else 'EMPTY'
            assert digest==r['inserted_md5'],(r['locus_id'],digest,r['inserted_md5'])
            r['canonical_fragment']=canonical_insert_key(long,start,length)
        summary[panel]=dict(input_loci=46,resolved=len(rows),missing=46-len(rows),
            classes=dict(collections.Counter(r['event_class'] for r in rows)),
            placements=dict(collections.Counter(r['placement_status'] for r in rows)))
    ea={r['event_id'].split('.')[-1]:r for r in read(root/'r6_first/armB_events.tsv')}
    eb={r['event_id'].split('.')[-1]:r for r in read(root/'tigr4_first/armB_events.tsv')}
    # Pair only after checking both exported sequence and biological carrier sets.
    comparison=[]
    for key,a in ea.items():
        b=eb[key]
        assert (a['insert_md5'],a['filled_genomes'],a['empty_genomes'])==(b['insert_md5'],b['filled_genomes'],b['empty_genomes'])
        x=results['r6_first'].get(key);y=results['tigr4_first'].get(key)
        row=dict(candidate=key,r6_resolved=bool(x),tigr4_resolved=bool(y),
                 same_class='.',same_insert_length='.',same_insert_md5='.',same_equivalent_fragment='.',
                 same_junction_frame='.',junction_L_shift='.',junction_R_shift='.',
                 r6_class=x['event_class'] if x else '.',tigr4_class=y['event_class'] if y else '.',
                 r6_insert_length=x['inserted_len'] if x else '.',tigr4_insert_length=y['inserted_len'] if y else '.',
                 r6_placement=x['placement_status'] if x else '.',tigr4_placement=y['placement_status'] if y else '.')
        if x and y:
            row.update(same_class=x['event_class']==y['event_class'],
                       same_insert_length=x['inserted_len']==y['inserted_len'],
                       same_insert_md5=x['inserted_md5']==y['inserted_md5'],
                       same_equivalent_fragment=x['canonical_fragment']==y['canonical_fragment'])
            frame=('pre_event_genome','pre_event_contig','pre_event_strand')
            same=all(x[k]==y[k] for k in frame)
            row['same_junction_frame']=same
            if same:
                row['junction_L_shift']=int(y['junction_L_abs'])-int(x['junction_L_abs'])
                row['junction_R_shift']=int(y['junction_R_abs'])-int(x['junction_R_abs'])
        comparison.append(row)
    common=[r for r in comparison if r['r6_resolved'] and r['tigr4_resolved']]
    summary['comparison']=dict(common_resolved=len(common),genome_allele_checks=sequences_checked,
        same_class=sum(r['same_class'] for r in common),same_insert_length=sum(r['same_insert_length'] for r in common),
        same_insert_md5=sum(r['same_insert_md5'] for r in common),
        same_equivalent_fragment=sum(r['same_equivalent_fragment'] for r in common),
        exact_same_junctions=sum(r['same_junction_frame'] and r['junction_L_shift']==0 and r['junction_R_shift']==0 for r in common))
    simple=[r for r in common if r['r6_class'] in ('insertion_target_retained','replacement','pure_insertion')]
    summary['shared_simple_representations']=dict(count=len(simple),
        same_fragment=sum(r['same_insert_md5'] for r in simple),
        same_junctions=sum(r['same_junction_frame'] and r['junction_L_shift']==0 and r['junction_R_shift']==0 for r in simple),
        unique_both=sum(r['r6_placement']=='unique' and r['tigr4_placement']=='unique' for r in simple))
    with (root/'layer2_comparison.tsv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(comparison[0]),delimiter='\t');w.writeheader();w.writerows(comparison)
    (root/'layer2_comparison.json').write_text(json.dumps(summary,indent=2))
    cols=list(comparison[0]);body=''.join('<tr>'+''.join('<td>'+html.escape(str(r[k]))+'</td>' for k in cols)+'</tr>' for r in comparison)
    (root/'layer2_comparison.html').write_text('<!doctype html><meta charset="utf-8"><title>Arm B to Layer 2</title><style>body{font:15px Arial;margin:30px}td,th{padding:9px;border-bottom:1px solid #ddd}table{border-collapse:collapse}</style><h1>Arm B → Layer 2: backbone-order comparison</h1><pre>'+html.escape(json.dumps(summary,indent=2))+'</pre><p>46 candidates in each order; unmodified Layer 2; 200 bp seed padding and 500 bp anchors. Missing loci are unresolved, not negatives. Class names describe sequence representations, not evolutionary direction. Padded allele lengths are not target motif lengths. Equivalent fragment comparison allows reverse complement and junction sliding within repeats. Sequence-export checks establish consistency, not unique locus identity or mechanistic precision.</p><p><a href="layer2_comparison.tsv">All 46 candidate comparisons</a></p><table><tr>'+''.join('<th>'+k+'</th>' for k in cols)+'</tr>'+body+'</table>',encoding='utf-8')
    print(json.dumps(summary,indent=2))
    print('Four previously shifted loci:')
    for r in comparison:
        if r['candidate'] in ('B000002','B000012','B000019','B000036'):print(json.dumps(r))


if __name__=='__main__':main()
