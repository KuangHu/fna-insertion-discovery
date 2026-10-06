import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from compare_pangraph_mummer import extract, distance


class PangraphComparisonTests(unittest.TestCase):
    def graph(self):
        return dict(paths={'0':dict(name='NC_003098.1',nodes=[1,2,3]),
                           '1':dict(name='NC_003028.3',nodes=[4,5])},
            nodes={str(i):dict(block_id=b,strand='+',position=p) for i,b,p in
                   [(1,10,[0,500]),(2,20,[500,1400]),(3,30,[1400,2000]),
                    (4,10,[0,500]),(5,30,[500,1100])]},
            blocks={'10':dict(alignments={})})

    def test_empty_versus_long_coordinates(self):
        rows,_=extract(self.graph())
        self.assertEqual(len(rows),1)
        self.assertEqual((rows[0]['r6_bp'],rows[0]['tigr4_bp'],rows[0]['signed_delta']),(900,0,-900))
        self.assertEqual(rows[0]['category'],'primary')

    def test_reoriented_anchor_pair_excluded(self):
        g=self.graph();g['nodes']['5']['strand']='-'
        rows,stats=extract(g)
        self.assertEqual(rows,[])
        self.assertEqual(stats['excluded_order_or_wrap'],1)

    def test_distant_second_genome_not_context_match(self):
        self.assertEqual(distance(10,20,15,25),0)
        self.assertGreater(distance(10,20,500,510),200)


if __name__=='__main__':unittest.main()
