"""Gate refactoring against the previously audited real-genome outputs."""
import csv
import json
from pathlib import Path
from local_debug import ROOT


def read(p):
    with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))


def main():
    checks={}
    for panel in ('r6_first','tigr4_first'):
        old=ROOT/'local_debug/armB'/panel
        new=ROOT/'local_debug/refactored'/panel
        state=json.loads((new/'summary.json').read_text())
        assert state['primary_status']==state['supplementary_status']=='complete'
        for filename in ('armB_events.tsv','armB_large_events.tsv','armB_pairs.tsv'):
            assert read(old/filename)==read(new/'armB'/filename),(panel,filename)
        for filename in ('armB_inserts.fna','armB_large_inserts.fna'):
            assert (old/filename).read_bytes()==(new/'armB'/filename).read_bytes(),(panel,filename)
        for filename in ('loci.tsv','alleles.tsv'):
            assert read(old/'layer2'/filename)==read(new/'layer2'/filename),(panel,filename)
        assert (old/'layer2/alleles.fna').read_bytes()==(new/'layer2/alleles.fna').read_bytes()
        assert len(read(new/'candidates.tsv'))==46
        assert sum(state['layer2_status_counts'].values())==46
        old_graph=json.loads((ROOT/'local_debug/pangraph_mummer'/(panel+'.json')).read_text())
        new_graph=json.loads((new/'supplementary/pangraph.json').read_text())
        assert old_graph==new_graph,(panel,'Pangraph graph')
        if panel=='r6_first':
            assert (ROOT/'local_debug/pangraph_mummer/all.coords').read_bytes()==(new/'supplementary/all.coords').read_bytes()
        assert any((new/'layer2_tmp').glob('allrec.*/*.paf'))
        checks[panel]=dict(primary_and_layer2_equal_to_baseline=True,pangraph_equal_to_baseline=True,
                           all_candidates_retained=True,status_counts=state['layer2_status_counts'])
    out=ROOT/'local_debug/refactored/verification.json'
    out.write_text(json.dumps(checks,indent=2));print(json.dumps(checks,indent=2))


if __name__=='__main__':main()
