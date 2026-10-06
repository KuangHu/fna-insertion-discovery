"""Independent assembly-level filled/short-allele audit of local Arm A calls.

Requires both inward anchor ends aligned, collinear same-contig placements,
and exactly one qualifying anchor pair. Missing mappings are never absences.
No evolutionary polarity or experimentally validated insertion is inferred.
"""
import collections
import csv
import html
import json
from pathlib import Path
import subprocess
from local_debug import ROOT, SOURCES, fasta


def read(p):
    with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))


def rc(s):return s.translate(str.maketrans('ACGTacgtNn','TGCAtgcaNn'))[::-1]


def main():
    out=ROOT/'local_debug/filled_empty';out.mkdir(parents=True,exist_ok=True)
    genomes={name:fasta(ROOT/'local_debug/genomes'/sp/(name+'.fna')) for name,sp,_,_ in SOURCES}
    species={name:sp for name,sp,_,_ in SOURCES}
    records=[];anchor=500
    with (out/'anchors.fna').open('w') as q,(out/'filled_alleles.fna').open('w') as filled:
        for name,sp,_,_ in SOURCES:
            for c in read(ROOT/'local_debug/armA'/(name+'_copies.tsv')):
                s,e=int(c['start']),int(c['end']);seq=genomes[name][c['contig']]
                if s<anchor or e+anchor>len(seq):continue
                ident='V%04d'%len(records)
                L,I,R=seq[s-anchor:s],seq[s:e],seq[e:e+anchor]
                assert L+I+R==seq[s-anchor:e+anchor]
                records.append(dict(id=ident,copy=c,source=name,L=L,I=I,R=R))
                q.write('>'+ident+'#L\n'+L+'\n>'+ident+'#R\n'+R+'\n')
                filled.write('>'+ident+' source='+name+' contig='+c['contig']+' interval='+str(s-anchor)+':'+str(e+anchor)+' frame=genomic\n'+L+I+R+'\n')
    hits={}
    for name,sp,_,_ in SOURCES:
        paf=out/(name+'.paf')
        with paf.open('w') as f,(out/(name+'.log')).open('w') as log:
            subprocess.run([str(ROOT/'env/bin/minimap2'),'-c','-x','asm20','-k','15','-w','5',
                '-N','100','-p','0.1','--secondary=yes','-t','2',
                str(ROOT/'local_debug/genomes'/sp/(name+'.fna')),str(out/'anchors.fna')],
                stdout=f,stderr=log,check=True)
        by=collections.defaultdict(list)
        for line in paf.read_text().splitlines():
            f=line.split('\t');ql,qs,qe=map(int,f[1:4]);side=f[0][-1]
            identity=int(f[9])/int(f[10]);cov=(qe-qs)/ql
            if identity<.95 or cov<.90:continue
            if (side=='L' and qe!=ql) or (side=='R' and qs!=0):continue
            by[f[0]].append(dict(contig=f[5],start=int(f[7]),end=int(f[8]),strand=f[4],identity=identity,cov=cov))
        hits[name]=by
    result=[]
    with (out/'comparison_alleles.fna').open('w') as alleles:
        for rec in records:
            ident,c,source=rec['id'],rec['copy'],rec['source']
            targets=[n for n in genomes if n!=source and species[n]==species[source]]
            if not targets:targets=[None]
            for target in targets:
                row=dict(id=ident,copy_id=c['copy_id'],family_id=c['family_id'],source_genome=source,
                    source_contig=c['contig'],source_start=c['start'],source_end=c['end'],
                    source_insert_bp=len(rec['I']),boundary_refined=c['boundary_refined'],
                    target_genome=target or '.',status='no_comparator',qualifying_pairs=0,
                    target_contig='.',target_start='.',target_end='.',target_strand='.',
                    observed_gap_bp='.',left_identity='.',right_identity='.',short_interval_sequence='.')
                if target:
                    pairs={}
                    for L in hits[target].get(ident+'#L',[]):
                        for R in hits[target].get(ident+'#R',[]):
                            if L['contig']!=R['contig'] or L['strand']!=R['strand']:continue
                            a,b=(L['end'],R['start']) if L['strand']=='+' else (R['end'],L['start'])
                            if not -50<=b-a<=10000:continue
                            pairs[L['contig'],a,b,L['strand']]=(L,R)
                    row['qualifying_pairs']=len(pairs)
                    row['status']='unresolved_anchors' if not pairs else 'ambiguous_placement'
                    if len(pairs)==1:
                        (ctg,a,b,strand),(L,R)=next(iter(pairs.items()));gap=b-a
                        row.update(target_contig=ctg,target_start=a,target_end=b,target_strand=strand,
                            observed_gap_bp=gap,left_identity=round(L['identity'],5),right_identity=round(R['identity'],5))
                        row['status']=('overlapping_anchors' if gap<0 else 'observed_zero_interval' if gap==0
                            else 'observed_short_interval' if gap<=50 and gap<=.2*len(rec['I'])
                            else 'similar_length_interval' if .8*len(rec['I'])<=gap<=1.2*len(rec['I'])
                            else 'other_length_interval')
                        if gap>=0:
                            seq=genomes[target][ctg];middle=seq[a:b]
                            lo,hi=min(L['start'],R['start']),max(L['end'],R['end'])
                            full=seq[lo:hi]
                            if strand=='-':middle,full=rc(middle),rc(full)
                            if gap<=50:row['short_interval_sequence']=middle or '(zero-length)'
                            alleles.write('>'+ident+' target='+target+' contig='+ctg+' interval='+str(lo)+':'+str(hi)+' strand='+strand+' status='+row['status']+'\n'+full+'\n')
                result.append(row)
    with (out/'verification.tsv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(result[0]),delimiter='\t');w.writeheader();w.writerows(result)
    summary=dict(source_copies_checked=len(records),counts=dict(collections.Counter(r['status'] for r in result)),
        method='500bp anchors; identity>=0.95; coverage>=0.90; both inward ends unclipped; unique collinear pair',
        limits='Assembly evidence only. Short intervals are preserved. No polarity, independent-event or RNA-guidance claims.')
    (out/'summary.json').write_text(json.dumps(summary,indent=2))
    cols=['copy_id','source_insert_bp','target_genome','status','observed_gap_bp','qualifying_pairs','boundary_refined']
    body=''.join('<tr>'+''.join('<td>'+html.escape(str(r[k]))+'</td>' for k in cols)+'</tr>' for r in sorted(result,key=lambda r:(r['status'],r['copy_id'])))
    (out/'index.html').write_text('<!doctype html><meta charset="utf-8"><title>Filled / short allele audit</title><style>body{font:16px Arial;margin:36px;color:#17354a}table{border-collapse:collapse}td,th{padding:9px;border-bottom:1px solid #ddd;text-align:left}th{background:#edf3f5}pre{background:#edf3f5;padding:16px}</style><h1>Filled / short allele audit</h1><p>Source: [500 bp host flank] [candidate insertion] [500 bp host flank]</p><p>Comparison: two homologous anchors on one continuous contig; the intervening sequence is extracted verbatim.</p><pre>'+html.escape(json.dumps(summary,indent=2))+'</pre><p>Zero/short intervals support an assembly-level long/short allele pair. Missing mappings are unresolved, never empty. Similar length does not establish sequence identity. Rows are candidate copies, not independent events.</p><table><tr>'+''.join('<th>'+k+'</th>' for k in cols)+'</tr>'+body+'</table>')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
