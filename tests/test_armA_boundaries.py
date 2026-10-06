import importlib.util
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
spec=importlib.util.spec_from_file_location('armA',ROOT/'scripts/20_armA_multicopy.py')
arm=importlib.util.module_from_spec(spec);spec.loader.exec_module(arm)


def paf(cigar,qe,te,qs=0,ts=0,strand='+'):
    return SimpleNamespace(tags={'cg':cigar},qs=qs,qe=qe,ts=ts,te=te,strand=strand)


class BoundaryTests(unittest.TestCase):
    def test_match_and_endpoints(self):
        r=paf('100M',120,150,20,50)
        self.assertEqual([arm.project_target_boundary(r,x) for x in [50,75,150]],[20,45,120])

    def test_deletion_and_insertion_offsets(self):
        r=paf('216M125D241M11I767M',1235,1349)
        self.assertEqual(arm.project_target_boundary(r,1300),1186)

    def test_insertion_junction_is_ambiguous(self):
        r=paf('20M5I80M',105,100)
        self.assertIsNone(arm.project_target_boundary(r,20))
        self.assertEqual(arm.project_target_boundary(r,21),26)

    def test_deletion_interior_is_unresolved(self):
        r=paf('20M5D80M',100,105)
        self.assertIsNone(arm.project_target_boundary(r,22))
        self.assertEqual(arm.project_target_boundary(r,20),20)
        self.assertEqual(arm.project_target_boundary(r,25),20)

    def test_no_extrapolation(self):
        r=paf('100M',120,150,20,50)
        self.assertIsNone(arm.project_target_boundary(r,49))
        self.assertIsNone(arm.project_target_boundary(r,151))

    def test_bad_or_missing_cigar(self):
        for cg in ['', '10S90M', '99M', '0M100M']:
            self.assertIsNone(arm.project_target_boundary(paf(cg,100,100),50))

    def test_equal_and_mismatch_ops(self):
        self.assertEqual(arm.project_target_boundary(paf('20=5X75=',100,100),30),30)

    def test_reverse_rejected(self):
        self.assertIsNone(arm.project_target_boundary(paf('100M',100,100,strand='-'),30))

    def test_refinement_integration_and_unaligned_status(self):
        ctx=[('ref','A'*1500,100,1100),('copy','A'*1500,100,1100),('nohit','A'*1500,100,1100)]
        def fake_run(cmd,stdout,**kwargs):
            stdout.write('Q1\t1500\t100\t1090\t+\tREP\t1500\t100\t1100\t990\t1000\t60\tcg:Z:500M10D490M\n')
        with tempfile.TemporaryDirectory() as tmp,patch.object(arm,'run',fake_run):
            anchors,stats=arm.refine_boundaries(ctx,0,str(Path(tmp)/'f'),'minimap2',1)
        self.assertEqual(anchors[1],(100,1090))
        self.assertEqual(anchors[2],(100,1100))
        self.assertEqual(stats['refined_indices'],{0,1})

    def test_one_export_coordinate_source(self):
        copies=[('ctg',100,1100),('ctg',2000,3000)]
        self.assertEqual(arm.resolved_copy(copies,{0:(110,1090)},0),('ctg',110,1090))
        self.assertEqual(arm.resolved_copy(copies,{0:(110,1090)},1),('ctg',2000,3000))

    def test_ambiguous_end_keeps_whole_nominal_interval(self):
        ctx=[('ref','A'*1500,100,1100),('copy','A'*1500,90,1120)]
        def fake_run(cmd,stdout,**kwargs):
            stdout.write('Q1\t1500\t100\t1105\t+\tREP\t1500\t100\t1100\t1000\t1005\t60\tcg:Z:5I1000M\n')
        with tempfile.TemporaryDirectory() as tmp,patch.object(arm,'run',fake_run):
            anchors,stats=arm.refine_boundaries(ctx,0,str(Path(tmp)/'f'),'minimap2',1)
        self.assertEqual(anchors[1],(90,1120))
        self.assertNotIn(1,stats['refined_indices'])
        self.assertEqual(stats['projection_unresolved'],1)


if __name__=='__main__':unittest.main()
