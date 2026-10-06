"""Synthetic end-to-end tests, independent of local reference downloads."""
import contextlib
import csv
import io
import json
from pathlib import Path
import random
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from classify_allele_evidence import rc, run


class AlleleModuleTests(unittest.TestCase):
    def setUp(self):
        rng=random.Random(154)
        self.L,self.I,self.R=[''.join(rng.choices('ACGT',k=n)) for n in (500,800,500)]
        self.genomes={'source':{'chr':self.L+self.I+'ATTATATAT'+self.R},
                      'target':{'chr':self.L+'ATTATATAT'+self.R}}
        self.row=dict(id='test',copy_id='copy1',source_genome='source',source_contig='chr',
            source_start='500',source_end='1309',source_insert_bp='809',boundary_refined='0',
            target_genome='target',target_contig='chr',target_start='500',target_end='509',
            target_strand='+',observed_gap_bp='9',qualifying_pairs='1',status='observed_short_interval')

    def execute(self,rows=None):
        with tempfile.TemporaryDirectory() as directory:
            with contextlib.redirect_stdout(io.StringIO()):
                summary=run(rows or [self.row],self.genomes,directory)
            with (Path(directory)/'uniform_evidence.tsv').open() as f:
                result=list(csv.DictReader(f,delimiter='\t'))
            seq=(Path(directory)/'uniform_middle_sequences.fna').read_text()
            return summary,result,seq

    def test_short_sequence_preserved_and_boundary_separate(self):
        _,rows,seq=self.execute()
        self.assertEqual(rows[0]['comparison_interval_sequence'],'ATTATATAT')
        self.assertIn('ATTATATAT',seq)
        self.assertEqual(rows[0]['evidence_status'],'supported_long_short')
        self.assertEqual(rows[0]['boundary_status'],'nominal_boundary')

    def test_reverse_genome_has_same_evidence_and_oriented_middle(self):
        _,forward,_=self.execute()
        self.genomes['target']['chr']=rc(self.genomes['target']['chr'])
        self.row['target_strand']='-'
        _,reverse,_=self.execute()
        for field in ('evidence_status','comparison_interval_sequence','dominant_gap_bp','gap_call_coverage'):
            self.assertEqual(forward[0][field],reverse[0][field])

    def test_anchor_overlap_kept_undefined_without_rejection(self):
        self.genomes['target']['chr']=self.L+self.R
        self.row.update(target_start='501',target_end='500',observed_gap_bp='-1',status='overlapping_anchors')
        _,rows,_=self.execute()
        self.assertEqual(rows[0]['comparison_interval_bp'],'.')
        self.assertEqual(rows[0]['boundary_status'],'anchor_overlap')
        self.assertEqual(rows[0]['evidence_status'],'supported_long_short')

    def test_unresolved_ambiguous_and_no_comparator_retained(self):
        rows=[dict(self.row,id='missing',qualifying_pairs='0',status='unresolved_anchors'),
              dict(self.row,id='ambiguous',qualifying_pairs='2',status='ambiguous_placement'),
              dict(self.row,id='alone',qualifying_pairs='0',status='no_comparator')]
        summary,result,_=self.execute(rows)
        self.assertEqual(summary['total'],3)
        self.assertEqual([r['evidence_status'] for r in result],['unresolved','unresolved','no_comparator'])

    def test_zero_middle_recorded_explicitly(self):
        self.genomes['target']['chr']=self.L+self.R
        self.row.update(target_end='500',observed_gap_bp='0',status='observed_zero_interval')
        _,rows,_=self.execute()
        self.assertEqual(rows[0]['comparison_interval_bp'],'0')
        self.assertEqual(rows[0]['comparison_interval_sequence'],'(zero-length)')

    def test_split_gap_is_not_reported_as_absence(self):
        self.genomes['target']['chr']=self.L+'ATTATATATAT'+self.R
        self.row.update(target_end='511',observed_gap_bp='11')
        _,rows,_=self.execute()
        self.assertEqual(rows[0]['evidence_status'],'no_dominant_long_short_support')
        self.assertGreater(int(rows[0]['dominant_gap_bp']),0)
        self.assertEqual(rows[0]['comparison_interval_sequence'],'ATTATATATAT')

    def test_invalid_inputs_fail_before_writing(self):
        variants=[[],[self.row,self.row], [dict(self.row,source_end='99999')],
                  [dict(self.row,target_strand='?')],[dict(self.row,observed_gap_bp='0')],
                  [dict(self.row,id='../escape')]]
        for rows in variants:
            with self.subTest(rows=rows),tempfile.TemporaryDirectory() as directory:
                out=Path(directory)/'out'
                with self.assertRaises(ValueError):run(rows,self.genomes,out)
                self.assertFalse(out.exists())

    def test_custom_manifest_cli(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            for name,contigs in self.genomes.items():
                (root/(name+'.fna')).write_text('>chr\n'+contigs['chr']+'\n')
            (root/'manifest.tsv').write_text('genome_id\tfasta_path\nsource\tsource.fna\ntarget\ttarget.fna\n')
            with (root/'calls.tsv').open('w',newline='') as f:
                writer=csv.DictWriter(f,fieldnames=list(self.row),delimiter='\t')
                writer.writeheader();writer.writerow(self.row)
            script=Path(__file__).resolve().parents[1]/'tools/classify_allele_evidence.py'
            subprocess.run([sys.executable,str(script),'--manifest',str(root/'manifest.tsv'),
                '--verification',str(root/'calls.tsv'),'--outdir',str(root/'result')],
                check=True,capture_output=True,text=True)
            summary=json.loads((root/'result/uniform_summary.json').read_text())
            self.assertEqual(summary['counts'],{'supported_long_short':1})


if __name__=='__main__':unittest.main()
