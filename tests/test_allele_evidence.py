import random
import sys
import unittest
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from classify_allele_evidence import compare, rc


class AlleleEvidenceTests(unittest.TestCase):
    def setUp(self):
        rng=random.Random(718)
        self.left,self.insert,self.right=[''.join(rng.choices('ACGT',k=n)) for n in (300,800,300)]

    def test_retained_short_sequence_supported(self):
        x=self.left+self.insert+'TTTTT'+self.right
        y=self.left+'TTTTT'+self.right
        result,_=compare(x,y,300,1105)
        self.assertEqual(result['evidence_status'],'supported_long_short')
        self.assertEqual(result['dominant_gap_bp'],800)

    def test_same_allele_not_supported(self):
        x=self.left+self.insert+self.right
        result,_=compare(x,x,300,1100)
        self.assertEqual(result['dominant_gap_bp'],0)
        self.assertEqual(result['evidence_status'],'no_dominant_long_short_support')

    def test_gap_outside_call_not_support(self):
        x=self.left+self.insert+self.right
        result,_=compare(x,self.left+self.right,0,250)
        self.assertEqual(result['evidence_status'],'no_dominant_long_short_support')

    def test_reverse_comparison_oriented_before_alignment(self):
        x=self.left+self.insert+self.right
        genomic=rc(self.left+self.right)
        result,_=compare(x,rc(genomic),300,1100)
        self.assertEqual(result['evidence_status'],'supported_long_short')


if __name__=='__main__':unittest.main()
