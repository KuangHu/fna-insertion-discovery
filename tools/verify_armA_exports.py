"""Check real-genome exports against source slices, coordinates and MD5."""
import csv
import hashlib
import json
from pathlib import Path
from local_debug import fasta, ROOT, SOURCES


def read(path):
    with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))


def rc(seq):return seq.translate(str.maketrans('ACGTacgtNn','TGCAtgcaNn'))[::-1]


def main():
    summary=[]
    for name,species,_,_ in SOURCES:
        prefix=ROOT/'local_debug/armA'/name
        fam=read(Path(str(prefix)+'_families.tsv'))
        copies=read(Path(str(prefix)+'_copies.tsv'))
        elements=fasta(Path(str(prefix)+'_elements.fna'))
        sites=fasta(Path(str(prefix)+'_sites.fna'))
        genome=fasta(ROOT/'local_debug/genomes'/species/(name+'.fna'))
        for c in copies:
            s,e=int(c['start']),int(c['end']);seq=genome[c['contig']]
            assert e-s==int(c['length'])
            left,element,right=seq[max(0,s-500):s],seq[s:e],seq[e:e+500]
            if c['strand']=='-':left,element,right=rc(right),rc(element),rc(left)
            assert sites[c['copy_id']]==left+element+right,c['copy_id']
            assert c['left_flank_200']==left[-200:]
            assert c['right_flank_200']==right[:200]
        for f in fam:
            ctg,bounds=f['rep_copy'].rsplit(':',1);s,e=map(int,bounds.split('-'))
            matches=[c for c in copies if c['family_id']==f['family_id'] and c['contig']==ctg and int(c['start'])==s and int(c['end'])==e]
            assert len(matches)==1,f['family_id']
            seq=genome[ctg][s:e]
            if matches[0]['strand']=='-':seq=rc(seq)
            assert elements[f['family_id']]==seq,f['family_id']
            assert hashlib.md5(seq.upper().encode()).hexdigest()==f['elem_md5']
        summary.append(dict(genome=name,representatives_checked=len(fam),copies_checked=len(copies),
            reverse_copies_checked=sum(c['strand']=='-' for c in copies),
            unrefined_copies=sum(c['boundary_refined']=='0' for c in copies),status='pass'))
    (ROOT/'local_debug/export_verification.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
