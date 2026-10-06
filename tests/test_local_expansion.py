"""CLI regression for information loss before record-level auditing."""
import csv
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]


@unittest.skipIf(os.name=='nt','Mock aligner uses a POSIX executable; run in WSL')
class TestExpansion(unittest.TestCase):
    def test_same_bin_hits_and_orientation_survive(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);(p/'genomes').mkdir()
            (p/'genomes/g.fna').write_text('>c\n'+'A'*3000+'\n')
            (p/'inserts.fna').write_text('>e\n'+'A'*100+'\n')
            (p/'cds.tsv').write_text('db_id\tcds_cluster_id\ne\tC1\n')
            (p/'db.tsv').write_text('db_id\tspecies\ne\tsource\n')
            paf=[]
            for start,strand in [(1000,'+'),(1020,'-')]:
                paf.append('\t'.join(map(str,['C1|e|source',100,0,100,strand,'g|c',3000,start,start+100,100,100,60])))
            fake=p/'aligner'
            fake.write_text('#!'+sys.executable+'\nprint('+repr('\n'.join(paf))+')\n')
            fake.chmod(0o755)
            cmd=[sys.executable,str(ROOT/'scripts/181_cross_species_expansion.py'),
                '--inserts',str(p/'inserts.fna'),'--cds',str(p/'cds.tsv'),'--db',str(p/'db.tsv'),
                '--target-dir',str(p/'genomes'),'--target-species','target',
                '--out',str(p/'result'),'--minimap2',str(fake)]
            subprocess.run(cmd,check=True,capture_output=True,text=True)
            with (p/'result_hits.tsv').open() as f:
                rows=list(csv.DictReader(f,delimiter='\t'))
            self.assertEqual(len(rows),2)
            self.assertEqual([r['strand'] for r in rows],['+','-'])
            self.assertEqual([r['query_start'] for r in rows],['0','0'])
            fake.write_text('#!'+sys.executable+'\nraise SystemExit(7)\n')
            failed=subprocess.run(cmd,capture_output=True,text=True)
            self.assertNotEqual(failed.returncode,0)
            self.assertIn('minimap2 failed',failed.stderr)


if __name__=='__main__':unittest.main()
