"""Conservative two-genome Pangraph path comparison and independent MUMmer audit."""
import csv
import collections
import html
import json
from local_debug import ROOT
from discovery_evidence import shared_intervals


def read(p):
    with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))


def table(p,rows):
    if not rows:return
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)


def extract(g):
    rows,stats=shared_intervals(g,['NC_003098.1','NC_003028.3'])
    translated=[dict(id=r['id'],r6_start=r['reference_start'],r6_end=r['reference_end'],
        tigr4_start=r['query_start'],tigr4_end=r['query_end'],r6_bp=r['reference_bp'],tigr4_bp=r['query_bp'],
        signed_delta=r['signed_delta'],category=r['category'],left_block=r['left_block'],right_block=r['right_block']) for r in rows]
    return translated,stats


def distance(a,b,c,d):
    return max(0,c-b,a-d)


def main():
    out=ROOT/'local_debug/pangraph_mummer';base=ROOT/'local_debug/armB'
    summary={};panels={}
    for panel in ('r6_first','tigr4_first'):
        rows,stats=extract(json.loads((out/(panel+'.json')).read_text()))
        table(out/(panel+'_intervals.tsv'),rows);panels[panel]=rows
        stats['interval_classes']=dict(collections.Counter(r['category'] for r in rows));summary[panel]=stats
    def key(r):return tuple(r[k] for k in ('r6_start','r6_end','tigr4_start','tigr4_end'))
    summary['order_identical_intervals']=set(map(key,panels['r6_first']))==set(map(key,panels['tigr4_first']))
    events=[r for r in read(base/'audit.tsv') if r['panel']=='r6_first' and r['category']=='primary']
    links=[];mapping={}
    for r in events:
        if r['filled_genome']=='spneumoniae_R6':
            s,e,u,v=map(int,(r['filled_start'],r['filled_end'],r['short_start'],r['short_end']))
        else:s,e,u,v=map(int,(r['short_start'],r['short_end'],r['filled_start'],r['filled_end']))
        mapping[r['event_id']]=(s,e,u,v)
        for p in panels['r6_first']:
            if p['category']=='below_500':continue
            if distance(s,e,p['r6_start'],p['r6_end'])<=200 and distance(u,v,p['tigr4_start'],p['tigr4_end'])<=200:
                links.append(dict(armB_event=r['event_id'],pangraph_interval=p['id'],pangraph_category=p['category'],
                    graph_delta=r['graph_length_difference'],pangraph_delta=abs(p['signed_delta']),
                    match_type='both_genome_context_within_200bp'))
    table(out/'candidate_links.tsv',links)
    linked={r['pangraph_interval'] for r in links};arm={r['armB_event'] for r in links}
    summary['context_comparison']=dict(armB_primary_with_context_match=len(arm),armB_primary_total=len(events),
        pangraph_primary_without_context_match=sum(r['category']=='primary' and r['id'] not in linked for r in panels['r6_first']),
        rule='Intervals overlap or are <=200 bp apart in BOTH genomes. Context match is NOT proof of identical events; nested/merged intervals possible.')
    coords=[]
    for line in (out/'all.coords').read_text().splitlines():
        f=line.split('\t')
        coords.append(dict(rs=int(f[0])-1,re=int(f[1]),qs=int(f[2])-1,qe=int(f[3]),
                           rlen=int(f[4]),qlen=int(f[5]),identity=float(f[6])))
    disputes=[]
    for suffix in ('B000002','B000012','B000019','B000036'):
        eid='r6_first.'+suffix;s,e,u,v=mapping[eid]
        region=[c for c in coords if distance(s,e,c['rs'],c['re'])<=2000]
        table(out/(suffix+'_all_regional_alignments.tsv'),region)
        good=[c for c in region if c['qs']<c['qe'] and min(c['rlen'],c['qlen'])>=200 and c['identity']>=95]
        left=[c for c in good if abs(c['re']-s)<=2000 and abs(c['qe']-u)<=2000 and c['re']<=e and c['qe']<=v]
        right=[c for c in good if abs(c['rs']-e)<=2000 and abs(c['qs']-v)<=2000 and c['rs']>=s and c['qs']>=u]
        pairs=[]
        for L in left:
            for R in right:
                if L['re']<=R['rs'] and L['qe']<=R['qs']:
                    score=abs(L['re']-s)+abs(L['qe']-u)+abs(R['rs']-e)+abs(R['qs']-v)
                    pairs.append((score,L,R))
        row=dict(event=eid,regional_alignment_blocks=len(region),candidate_flank_pairs=len(pairs),
                 r6_gap_start='.',r6_gap_end='.',tigr4_gap_start='.',tigr4_gap_end='.',
                 independent_gap_delta='.',left_identity='.',right_identity='.',status='unresolved')
        if pairs:
            _,L,R=min(pairs,key=lambda z:z[0])
            row.update(r6_gap_start=L['re'],r6_gap_end=R['rs'],tigr4_gap_start=L['qe'],tigr4_gap_end=R['qs'],
                       independent_gap_delta=(R['qs']-L['qe'])-(R['rs']-L['re']),
                       left_identity=L['identity'],right_identity=R['identity'],
                       status='local_collinear_pair; boundaries_are_alignment_extents')
        disputes.append(row)
    table(out/'disputed_mummer.tsv',disputes)
    summary['mummer_disputes']=disputes
    (out/'summary.json').write_text(json.dumps(summary,indent=2))
    (out/'report.html').write_text('<!doctype html><meta charset="utf-8"><title>Pangraph and MUMmer comparison</title><style>body{font:16px Arial;margin:35px;max-width:1200px}pre{white-space:pre-wrap;background:#eef3f5;padding:20px}</style><h1>Pangraph + independent MUMmer alignment</h1><p>Pangraph 1.4.0, circular genomes, default alignment settings, reconstruction verification passed. MUMmer 3.23 nucmer --maxmatch; both raw and one-to-one filtered alignments retained.</p><p>Comparison uses adjacent single-copy shared blocks in matching orientation. Origin-crossing and rearranged anchor pairs are not classified. Within-block variants are not independently called. These limits prevent a whole-tool sensitivity claim.</p><pre>'+html.escape(json.dumps(summary,indent=2))+'</pre><p><a href="r6_first_intervals.tsv">All Pangraph intervals</a> · <a href="candidate_links.tsv">Context links to Arm B</a> · <a href="disputed_mummer.tsv">Four disputed loci</a></p><p>MUMmer flank pairs are selected near BOTH known genome loci using all alignments; they are independent alignment evidence, not independent discovery or proof of unique placement. Repeat hits can create misleading rearrangement labels in the one-to-one filter.</p>',encoding='utf-8')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
