#!/usr/bin/env python3
"""Run locked Layer 2 on both local Arm B orders, using documented 200 bp pad."""
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1]


def main():
    root=ROOT/'local_debug/armB'
    script=ROOT/'scripts/70_allele_reconstructor.py'
    digest=hashlib.sha256(script.read_bytes()).hexdigest()
    for panel in ('r6_first','tigr4_first'):
        with (root/panel/'armB_events.tsv').open() as f:
            events=list(csv.DictReader(f,delimiter='\t'))
        seeds=root/panel/'layer2_seeds.tsv'
        with seeds.open('w',newline='') as f:
            w=csv.writer(f,delimiter='\t');w.writerow(['locus_id','seed_genome','contig','start','end'])
            for r in events:
                left,right=int(r['junction_L']),int(r['junction_R'])
                if left<200 or right<left:raise ValueError('Unusable seed: '+r['event_id'])
                w.writerow([r['event_id'],r['empty_ref_genome'],r['empty_ref_contig'],left-200,right+200])
        cmd=[sys.executable,str(script),'--manifest',str(root/(panel+'.manifest')),
             '--loci',str(seeds),'--outdir',str(root/panel/'layer2'),'--threads','2',
             '--minimap2',str(ROOT/'env/bin/minimap2')]
        with (root/panel/'layer2.log').open('w') as log:
            subprocess.run(cmd,stdout=log,stderr=log,check=True)
        print(panel+' Layer 2 complete',flush=True)
    assert digest==hashlib.sha256(script.read_bytes()).hexdigest()
    (root/'layer2_run.json').write_text(json.dumps(dict(script_sha256=digest,pad_bp=200,
        layer2_modified=False,parameters='defaults; 2 threads; local minimap2'),indent=2))


if __name__=='__main__':main()
