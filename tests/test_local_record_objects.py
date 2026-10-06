import importlib.util
from pathlib import Path
import random
import unittest

spec=importlib.util.spec_from_file_location('objects',Path(__file__).resolve().parents[1]/'scripts/local_record_objects.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


def row(start=100,end=1100,seed='seed1',family='fam1',**kwargs):
    r=dict(genome='g',contig='c',elem_start=str(start),elem_end=str(end),
        source_event=seed,cds_cluster=family,pct_ident='.99',q_cov='.95',
        flank_complete='yes',left_flank='A'*60,right_flank='C'*60,strand='+')
    r.update(kwargs);return r


class TestObjects(unittest.TestCase):
    def test_shuffle_stable_even_score_ties(self):
        rows=[row(),row(102,1105,'seed2','fam2'),row(1500,2500,'seed3')]
        expected=m.summarize(rows)[0]
        for i in range(20):
            shuffled=rows[:];random.Random(i).shuffle(shuffled)
            self.assertEqual(m.summarize(shuffled)[0],expected)

    def test_nested_hit_does_not_split_surrounding_group(self):
        groups=m.group_hits([row(),row(101,300,'nested'),row(102,1102,'outer2')])
        self.assertEqual(sorted(map(len,groups)),[1,2])

    def test_overlap_chain_not_merged(self):
        groups=m.group_hits([row(0,1000),row(90,1090,'b'),row(180,1180,'c')])
        self.assertEqual(sorted(map(len,groups)),[1,2])

    def test_flank_prefix_collision_avoided(self):
        rows=[row(),row(2000,3000,'other',left_flank='A'*12+'G'*48)]
        out,_=m.summarize(rows)
        self.assertNotEqual(out[0]['exact_flank_key'],out[1]['exact_flank_key'])

    def test_incomplete_or_ambiguous_flank_not_keyed(self):
        self.assertIsNone(m.flank_key(row(flank_complete='no')))
        self.assertIsNone(m.flank_key(row(left_flank='N'*60)))

    def test_seed_count_not_family_count(self):
        out,groups=m.summarize([row(),row(seed='second')])
        self.assertEqual(out[0]['n_seeds'],2)
        self.assertEqual(out[0]['n_families'],1)
        self.assertEqual(sum(map(len,groups)),2)

    def test_query_orientation_preserved_without_splitting_locus(self):
        out,_=m.summarize([row(),row(seed='reverse',strand='-')])
        self.assertEqual(len(out),1)
        self.assertEqual(out[0]['strands'],'+,-')

    def test_invalid_interval(self):
        with self.assertRaises(ValueError):m.group_hits([row(10,10)])


if __name__=='__main__':unittest.main()
