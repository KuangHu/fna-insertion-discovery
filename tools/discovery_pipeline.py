#!/usr/bin/env python3
"""Primary discovery: Arm B -> locked Layer 2; optional supplementary engines."""
import argparse
import collections
import csv
import hashlib
import html
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from Bio import SeqIO
from discovery_evidence import candidate_row, padded_seed, shared_intervals

ROOT=Path(__file__).resolve().parents[1]


def read(path):
    with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))


def write_table(path,rows,fields=None):
    with path.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=fields or list(rows[0]),delimiter='\t')
        writer.writeheader();writer.writerows(rows)


def executable(value):
    path=shutil.which(value)
    if path is None:raise ValueError('Executable not found: '+value)
    return str(Path(path).resolve())


def inputs(manifest):
    paths=[(manifest.parent/line.strip()).resolve() for line in manifest.read_text().splitlines() if line.strip()]
    if len(paths)<2 or len(set(paths))!=len(paths):
        raise ValueError('Manifest requires at least two distinct assembly files')
    records={};lengths={}
    for p in paths:
        name=p.name
        for suffix in ('.fna.gz','.fa.gz','.fasta.gz','.fna','.fa','.fasta'):
            if name.endswith(suffix):name=name[:-len(suffix)];break
        if name in lengths:raise ValueError('Duplicate assembly basename: '+name)
        if p.suffix=='.gz':raise ValueError('Use uncompressed FASTA for this runner')
        with p.open() as fasta:
            rr=SeqIO.to_dict(SeqIO.parse(fasta,'fasta'))
        if not rr:raise ValueError('Empty assembly: '+str(p))
        records[name]=rr;lengths[name]={k:len(v) for k,v in rr.items()}
    return paths,records,lengths


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--manifest',type=Path,required=True,help='One uncompressed FASTA path per line; relative to manifest')
    ap.add_argument('--outdir',type=Path,required=True,help='New or empty output directory')
    ap.add_argument('--panel-id',default='panel')
    ap.add_argument('--threads',type=int,default=2)
    for name in ('minigraph','gfatools','minimap2','samtools','pangraph','nucmer','delta-filter','show-coords','show-diff'):
        ap.add_argument('--'+name,default=name)
    ap.add_argument('--supplement',action='store_true',help='Pangraph + MUMmer evidence; requires exactly two single-contig assemblies')
    ap.add_argument('--circular',action='store_true',help='Declare circular genomes to supplementary Pangraph')
    args=ap.parse_args()
    try:
        paths,records,lengths=inputs(args.manifest.resolve())
        if args.threads<1:raise ValueError('threads must be positive')
        if not args.panel_id or any(c not in 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-' for c in args.panel_id):
            raise ValueError('panel-id must contain only letters, digits, underscore or hyphen')
        required=['minigraph','gfatools','minimap2','samtools']
        if args.supplement:
            if len(paths)!=2 or any(len(r)!=1 for r in records.values()):
                raise ValueError('Supplementary comparison currently requires two single-contig assemblies')
            names=[next(iter(r)) for r in records.values()]
            if len(set(names))!=2:raise ValueError('Pangraph requires distinct FASTA record IDs across the two assemblies')
            required+=['pangraph','nucmer','delta_filter','show_coords','show_diff']
        bins={name:executable(getattr(args,name)) for name in required}
        out=args.outdir.resolve()
        if out.exists() and any(out.iterdir()):raise ValueError('Output directory must be empty; use a new directory')
    except (ValueError,OSError,KeyError) as error:ap.error(str(error))
    out.mkdir(parents=True,exist_ok=True)
    logs=out/'logs';logs.mkdir()
    manifest=out/'manifest.txt';manifest.write_text(''.join(str(p)+'\n' for p in paths))
    provenance=dict(input_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
                    executable_sha256={k:hashlib.sha256(Path(v).read_bytes()).hexdigest() for k,v in bins.items()},
                    scripts_sha256={name:hashlib.sha256((ROOT/'scripts'/name).read_bytes()).hexdigest()
                                    for name in ('30_armB_graph.py','70_allele_reconstructor.py')},commands=[])
    summary=dict(primary_status='running',supplementary_status='not_requested',panel_id=args.panel_id)
    def save():
        (out/'provenance.json').write_text(json.dumps(provenance,indent=2))
        (out/'summary.json').write_text(json.dumps(summary,indent=2))
    def run(label,cmd,stdout=None,env=None):
        provenance['commands'].append(dict(stage=label,argv=list(map(str,cmd))))
        save()
        with (logs/(label+'.log')).open('w') as log:
            if stdout:
                with stdout.open('w') as f:subprocess.run(cmd,stdout=f,stderr=log,env=env,check=True)
            else:subprocess.run(cmd,stdout=log,stderr=log,env=env,check=True)
        print(label+' complete',flush=True)
    try:
        environment=os.environ.copy()
        environment['PATH']=str(Path(bins['samtools']).parent)+os.pathsep+environment.get('PATH','')
        arm=out/'armB'
        run('armB',[sys.executable,str(ROOT/'scripts/30_armB_graph.py'),'--manifest',str(manifest),
            '--outdir',str(arm),'--cluster-id',args.panel_id,'--threads',str(args.threads),
            '--minigraph',bins['minigraph'],'--gfatools',bins['gfatools']],env=environment)
        events=read(arm/'armB_events.tsv');seeds=[];statuses={}
        for event in events:
            seed,status=padded_seed(event,lengths)
            statuses[event['event_id']]=status
            if seed:seeds.append(seed)
        seedfile=out/'layer2_seeds.tsv'
        write_table(seedfile,seeds,['locus_id','seed_genome','contig','start','end'])
        layer=out/'layer2'
        if seeds:
            temporary=out/'layer2_tmp';temporary.mkdir()
            environment['TMPDIR']=str(temporary)
            run('layer2',[sys.executable,str(ROOT/'scripts/70_allele_reconstructor.py'),'--manifest',str(manifest),
                '--loci',str(seedfile),'--outdir',str(layer),'--threads',str(args.threads),'--minimap2',bins['minimap2']],env=environment)
            reconstructed={r['locus_id']:r for r in read(layer/'loci.tsv')}
        else:reconstructed={}
        candidates=[candidate_row(r,statuses[r['event_id']],reconstructed.get(r['event_id'])) for r in events]
        if candidates:write_table(out/'candidates.tsv',candidates)
        else:write_table(out/'candidates.tsv',[],['event_id','layer2_status','seed_status'])
        summary.update(primary_status='complete',primary_candidates=len(events),
            large_events=len(read(arm/'armB_large_events.tsv')),seeded=len(seeds),
            layer2_status_counts=dict(collections.Counter(r['layer2_status'] for r in candidates)))
        save()
        if args.supplement:
            summary['supplementary_status']='running';save()
            supp=out/'supplementary';supp.mkdir()
            cmd=[bins['pangraph'],'build','--verify','-j',str(args.threads),'-o',str(supp/'pangraph.json')]
            if args.circular:cmd.append('--circular')
            run('pangraph',cmd+list(map(str,paths)))
            intervals,stats=shared_intervals(json.loads((supp/'pangraph.json').read_text()),names)
            if intervals:write_table(supp/'pangraph_intervals.tsv',intervals)
            else:write_table(supp/'pangraph_intervals.tsv',[],['id','category'])
            summary['pangraph']=dict(stats,interval_classes=dict(collections.Counter(r['category'] for r in intervals)))
            run('nucmer',[bins['nucmer'],'--maxmatch','-p',str(supp/'pair')]+list(map(str,paths)))
            run('delta_filter',[bins['delta_filter'],'-1',str(supp/'pair.delta')],supp/'pair.1delta')
            for label,delta in [('all','pair.delta'),('one_to_one','pair.1delta')]:
                run('coords_'+label,[bins['show_coords'],'-THrcl',str(supp/delta)],supp/(label+'.coords'))
                run('diff_'+label,[bins['show_diff'],'-H',str(supp/delta)],supp/(label+'.diff'))
            summary['supplementary_status']='complete'
        save()
        body=''.join('<tr>'+''.join('<td>'+html.escape(str(r[k]))+'</td>' for k in
            ('event_id','insert_size','insert_len_direct','layer2_status','layer2_event_class','layer2_inserted_len'))+'</tr>' for r in candidates)
        (out/'report.html').write_text('<!doctype html><meta charset="utf-8"><title>Discovery pipeline</title><style>body{font:16px Arial;margin:30px}td,th{padding:10px;border-bottom:1px solid #ddd}pre{white-space:pre-wrap}</style><h1>Arm B → Layer 2</h1><pre>'+html.escape(json.dumps(summary,indent=2))+'</pre><p>All primary candidates are retained, including unseedable, unresolved and ambiguous loci. Graph length difference, exported fragment length and reconstructed difference length are separate. Pangraph and MUMmer remain supplementary evidence; their output does not filter primary candidates. Large events are parked separately. Sequence representation does not establish evolutionary polarity or biochemical junction accuracy.</p><p><a href="candidates.tsv">Candidate table</a> · <a href="provenance.json">Run provenance</a></p><table><tr><th>Candidate</th><th>Graph delta</th><th>Export length</th><th>Layer 2 status</th><th>Representation</th><th>Reconstructed length</th></tr>'+body+'</table>',encoding='utf-8')
    except Exception as error:
        summary['error']=str(error)
        if summary['primary_status']!='complete':summary['primary_status']='failed'
        else:summary['supplementary_status']='failed'
        save();raise
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
