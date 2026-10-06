#!/usr/bin/env python3
"""Run real two-genome Arm B panels in both backbone orders under WSL."""
import json
import subprocess
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]


def main():
    out=ROOT/'local_debug/armB'
    out.mkdir(parents=True,exist_ok=True)
    genomes=[ROOT/'local_debug/genomes/spneumoniae'/('spneumoniae_'+n+'.fna') for n in ('R6','TIGR4')]
    versions={n:subprocess.check_output(['git','-C',str(ROOT/'env/src'/n),'rev-parse','HEAD'],text=True).strip()
              for n in ('minigraph','gfatools')}
    (out/'dependency_commits.json').write_text(json.dumps(versions,indent=2))
    for label,paths in [('r6_first',genomes),('tigr4_first',genomes[::-1])]:
        manifest=out/(label+'.manifest')
        manifest.write_text(''.join(str(p)+'\n' for p in paths))
        cmd=[sys.executable,str(ROOT/'scripts/30_armB_graph.py'),'--manifest',str(manifest),
             '--outdir',str(out/label),'--cluster-id',label,'--threads','2',
             '--minigraph',str(ROOT/'env/bin/minigraph'),'--gfatools',str(ROOT/'env/bin/gfatools')]
        with (out/(label+'.log')).open('w') as log:
            subprocess.run(cmd,stdout=log,stderr=log,check=True)
        print(label+' complete',flush=True)


if __name__=='__main__':main()
