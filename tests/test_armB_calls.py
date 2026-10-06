"""Protect the distinction between graph path length and genomic span."""
import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
spec=importlib.util.spec_from_file_location('armb_calls',ROOT/'scripts/30_armB_graph.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


class ArmBCallTests(unittest.TestCase):
    def parse(self,record):
        with tempfile.TemporaryDirectory() as directory:
            p=Path(directory)/'calls.bed'
            p.write_text('chr\t100\t200\t>s1\t>s2\t'+record+'\n')
            return module.parse_call(p)[('chr',100,200)]

    def test_real_forward_record_keeps_distinct_lengths(self):
        r=self.parse('>s3:1443:+:NC_003098.1:40181:41626')
        self.assertEqual(r[:5],(1443,'NC_003098.1',40181,41626,'+'))
        self.assertNotEqual(r[0],r[3]-r[2])

    def test_zero_graph_path_still_has_genomic_span(self):
        r=self.parse('*:0:+:NC_003098.1:34548:34556')
        self.assertEqual(r[0],0)
        self.assertEqual(r[3]-r[2],8)

    def test_reverse_orientation_retained(self):
        r=self.parse('<s3:100:-:chr2:900:1002')
        self.assertEqual(r[:5],(100,'chr2',900,1002,'-'))

    def test_missing_call_not_zero_length(self):
        self.assertIsNone(self.parse('*')[0])


if __name__=='__main__':unittest.main()
