"""Independent local discovery and whole-genome alignment on the same inputs."""
from pathlib import Path
import subprocess
import json

ROOT=Path(__file__).resolve().parents[1]


def main():
    out=ROOT/'local_debug/pangraph_mummer';out.mkdir(parents=True,exist_ok=True)
    paths=[ROOT/'local_debug/genomes/spneumoniae'/('spneumoniae_'+n+'.fna') for n in ('R6','TIGR4')]
    commands=[]
    def run(label,cmd):
        commands.append(cmd)
        with (out/(label+'.log')).open('w') as log:
            subprocess.run(cmd,stdout=log,stderr=log,check=True)
        print(label+' complete',flush=True)
    for label,seqs in [('r6_first',paths),('tigr4_first',paths[::-1])]:
        run('pangraph_'+label,[str(ROOT/'env/bin/pangraph'),'build','--verify','--circular','-j','2',
            '-o',str(out/(label+'.json'))]+list(map(str,seqs)))
    run('nucmer',['nucmer','--maxmatch','-p',str(out/'r6_tigr4')]+list(map(str,paths)))
    with (out/'r6_tigr4.1delta').open('w') as f:
        subprocess.run(['delta-filter','-1',str(out/'r6_tigr4.delta')],stdout=f,check=True)
    for label,delta in [('all','r6_tigr4.delta'),('one_to_one','r6_tigr4.1delta')]:
        with (out/(label+'.coords')).open('w') as f:
            subprocess.run(['show-coords','-THrcl',str(out/delta)],stdout=f,check=True)
        with (out/(label+'.diff')).open('w') as f:
            subprocess.run(['show-diff','-H',str(out/delta)],stdout=f,check=True)
    (out/'commands.json').write_text(json.dumps(commands,indent=2))


if __name__=='__main__':main()
