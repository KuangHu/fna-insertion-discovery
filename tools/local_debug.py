#!/usr/bin/env python3
"""Download three small public references and run a local, no-SLURM smoke test.
Run from WSL: python3 tools/local_debug.py
No claim of validated mobile elements or guide-dependent activity.
"""
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import urllib.request

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'local_debug'
SOURCES=[
 ('spneumoniae_R6','spneumoniae','GCF_000007045.1_ASM704v1','000/007/045'),
 ('spneumoniae_TIGR4','spneumoniae','GCF_000006885.1_ASM688v1','000/006/885'),
 ('ecoli_K12','ecoli','GCF_000005845.2_ASM584v2','000/005/845')]


def fasta(path):
    seqs={}; name=None
    for line in path.read_text().splitlines():
        if line.startswith('>'):
            name=line[1:].split()[0];seqs[name]=''
        elif name:
            seqs[name]+=line.strip()
    return seqs


def table(path,fields,rows):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)


def main():
    os.chdir(ROOT)
    for d in ['logs','armA','expansion','seeds']:
        (OUT/d).mkdir(parents=True,exist_ok=True)
    manifest=[]; runs=[]
    def run(label,args):
        start=time.monotonic()
        with (OUT/'logs'/(label+'.log')).open('w') as f:
            subprocess.run(args,stdout=f,stderr=subprocess.STDOUT,check=True)
        runs.append(dict(label=label,wall_seconds=round(time.monotonic()-start,3),command=args))
        print(label+' complete',flush=True)
    for name,species,assembly,shard in SOURCES:
        folder=OUT/'genomes'/species;folder.mkdir(parents=True,exist_ok=True)
        dest=folder/(name+'.fna')
        url='https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/'+shard+'/'+assembly+'/'+assembly+'_genomic.fna.gz'
        if not dest.exists():
            print('Downloading '+assembly,flush=True)
            with urllib.request.urlopen(url,timeout=120) as response:
                compressed=response.read()
            sequence=gzip.decompress(compressed)
            if not sequence.startswith(b'>'):
                raise ValueError('Not FASTA')
            dest.write_bytes(sequence)
        seqs=fasta(dest)
        manifest.append(dict(name=name,species=species,assembly=assembly,url=url,
            sha256=hashlib.sha256(dest.read_bytes()).hexdigest(),contigs=len(seqs),
            bases=sum(map(len,seqs.values()))))
        run('armA_'+name,[sys.executable,'scripts/20_armA_multicopy.py',str(dest),
            '--out',str(OUT/'armA'/name),'--threads','2','--minimap2',str(ROOT/'env/bin/minimap2')])
    (OUT/'download_manifest.json').write_text(json.dumps(manifest,indent=2))
    # Debug-only family labels, NOT production CDS clusters. Use all Arm A
    # representatives to exercise real hits without fetching the HPC catalogue.
    cds=[]; db=[]; seq_records=[]
    for name,species,_,_ in SOURCES:
        for old,seq in fasta(OUT/'armA'/(name+'_elements.fna')).items():
            ident=name+'.'+old
            cds.append(dict(db_id=ident,cds_cluster_id='DEBUG%04d'%len(cds)))
            db.append(dict(db_id=ident,species=species))
            seq_records.append('>'+ident+'\n'+seq+'\n')
    if not cds:
        raise RuntimeError('No Arm A examples found')
    (OUT/'seeds/inserts.fna').write_text(''.join(seq_records))
    table(OUT/'seeds/cds.tsv',['db_id','cds_cluster_id'],cds)
    table(OUT/'seeds/db.tsv',['db_id','species'],db)
    for species in ['spneumoniae','ecoli']:
        prefix=str(OUT/'expansion'/species)
        run('expand_'+species,[sys.executable,'scripts/181_cross_species_expansion.py',
            '--inserts',str(OUT/'seeds/inserts.fna'),'--cds',str(OUT/'seeds/cds.tsv'),
            '--db',str(OUT/'seeds/db.tsv'),'--target-dir',str(OUT/'genomes'/species),
            '--target-species',species,'--out',prefix,'--threads','2',
            '--minimap2',str(ROOT/'env/bin/minimap2')])
        run('objects_'+species,[sys.executable,'scripts/local_record_objects.py',
            '--hits',prefix+'_hits.tsv','--out',prefix+'_objects'])
    (OUT/'run_manifest.json').write_text(json.dumps(dict(debug_seed_count=len(cds),runs=runs),indent=2))
    # Every accepted hit must survive in the audit's companion table.
    for species in ['spneumoniae','ecoli']:
        prefix=OUT/'expansion'/species
        def rows(path):
            with open(path) as f:
                return list(csv.DictReader(f,delimiter='\t'))
        original=rows(str(prefix)+'_hits.tsv')
        retained=rows(str(prefix)+'_objects.hits.tsv')
        for r in retained:
            r.pop('site_id')
        assert sorted(map(lambda r:json.dumps(r,sort_keys=True),original)) == sorted(map(lambda r:json.dumps(r,sort_keys=True),retained))
    print('Finished: '+str(OUT),flush=True)


if __name__=='__main__':
    main()
