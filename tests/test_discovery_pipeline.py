"""Guard candidate retention and seed coordinates at the workflow boundary."""
import sys
import tempfile
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from discovery_evidence import candidate_row, padded_seed
from discovery_pipeline import inputs


class DiscoveryPipelineTests(unittest.TestCase):
    def event(self):
        return dict(event_id='p.B000001',empty_ref_genome='g',empty_ref_contig='c',
                    junction_L='700',junction_R='705')

    def test_retained_short_interval_is_preserved_by_padding(self):
        seed,status=padded_seed(self.event(),{'g':{'c':2000}})
        self.assertEqual(status,'seeded')
        self.assertEqual((seed['start'],seed['end']),(500,905))

    def test_contig_edge_becomes_explicit_unseedable_row(self):
        event=self.event();event['junction_L']='100'
        seed,status=padded_seed(event,{'g':{'c':2000}})
        self.assertIsNone(seed)
        row=candidate_row(event,status,None)
        self.assertEqual(row['event_id'],event['event_id'])
        self.assertEqual(row['layer2_status'],'unseedable')

    def test_missing_and_ambiguous_loci_remain_candidates(self):
        self.assertEqual(candidate_row(self.event(),'seeded',None)['layer2_status'],'unresolved')
        row=candidate_row(self.event(),'seeded',dict(placement_status='ambiguous',inserted_len='912'))
        self.assertEqual(row['layer2_status'],'ambiguous')
        self.assertEqual(row['layer2_inserted_len'],'912')

    def test_invalid_frame_cannot_be_seeded(self):
        event=self.event();event['junction_R']='699'
        self.assertEqual(padded_seed(event,{'g':{'c':2000}}),(None,'invalid_junction_frame'))

    def test_manifest_paths_resolve_relative_to_manifest(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            for n in ('a','b'):(root/(n+'.fna')).write_text('>contig\nACGT\n')
            (root/'manifest').write_text('a.fna\nb.fna\n')
            paths,records,lengths=inputs(root/'manifest')
            self.assertEqual(paths,[root/'a.fna',root/'b.fna'])
            self.assertEqual(lengths,{'a':{'contig':4},'b':{'contig':4}})

    def test_duplicate_assemblies_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);(root/'manifest').write_text('a.fna\na.fna\n')
            with self.assertRaises(ValueError):inputs(root/'manifest')


if __name__=='__main__':unittest.main()
