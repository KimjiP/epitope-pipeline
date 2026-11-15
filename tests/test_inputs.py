"""Tests for what run_pipeline.py reads from UniProt and the PDB: which structures it keeps, which ids it
accepts, and how it writes a download. No network: the two calls to the outside are replaced.

Run:  python3 -I -m unittest discover -s tests -v
"""
import json
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'pipeline'))

import run_pipeline as rp          # noqa: E402


def loc(a, b):
    return {'start': {'value': a}, 'end': {'value': b}}


# A receptor that crosses the membrane several times, the way CCR5 and CD20 do
FAKE = {
    'sequence': {'value': 'M' * 300},
    'proteinDescription': {'recommendedName': {'fullName': {'value': 'Fake receptor'}}},
    'uniProtKBCrossReferences': [
        {'database': 'PDB', 'id': '1ABC', 'properties': [{'key': 'Chains', 'value': 'A=20-60'}]},
        {'database': 'PDB', 'id': '../x', 'properties': [{'key': 'Chains', 'value': 'A=20-60'}]},
        {'database': 'PDB', 'id': '2XYZ', 'properties': [{'key': 'Chains', 'value': 'A/B=150-200'}]},
        {'database': 'PDB', 'id': '3DEF', 'properties': [{'key': 'Chains', 'value': 'A=100-130'}]},
        {'database': 'PDB', 'id': '4GHI', 'properties': []},
        {'database': 'Pfam', 'id': 'PF00001'}],
    'features': [
        {'type': 'Signal', 'location': loc(1, 19)},
        {'type': 'Topological domain', 'description': 'Extracellular', 'location': loc(20, 60)},
        {'type': 'Topological domain', 'description': 'Cytoplasmic', 'location': loc(61, 99)},
        {'type': 'Topological domain', 'description': 'Extracellular', 'location': loc(140, 210)},
        {'type': 'Repeat', 'description': 'Approximate', 'location': loc(75, 300)},
        {'type': 'Domain', 'description': 'Ig-like V-type', 'location': loc(30, 140)}],
}


class UniProt(unittest.TestCase):
    def setUp(self):
        self._get = rp.get
        rp.get = lambda url, timeout=60: json.dumps(FAKE).encode()

    def tearDown(self):
        rp.get = self._get

    def test_every_extracellular_segment_is_kept(self):
        seq, name, pdbs, doms, sig, regions = rp.uniprot_entry('P00000')
        self.assertEqual(regions, [(20, 60), (140, 210)])      # the first version kept only (140, 210)
        self.assertEqual(sig, 19)

    def test_only_real_pdb_ids_are_kept(self):
        pdbs = rp.uniprot_entry('P00000')[2]
        self.assertEqual([p['id'] for p in pdbs], ['1ABC', '2XYZ', '3DEF', '4GHI'])

    def test_a_domain_called_approximate_is_not_a_name(self):
        doms = rp.uniprot_entry('P00000')[3]
        self.assertEqual([d[0] for d in doms], ['Ig-like-V-type-1'])

    def test_structures_are_kept_when_they_touch_any_segment(self):
        seq, name, pdbs, doms, sig, regions = rp.uniprot_entry('P00000')
        keep = [p['id'] for p in rp.in_regions(pdbs, regions)]
        self.assertEqual(keep, ['1ABC', '2XYZ'])               # 3DEF (100-130) is in the cytoplasmic gap; 4GHI has no range
        last_only = [p['id'] for p in rp.in_regions(pdbs, [regions[-1]])]
        self.assertEqual(last_only, ['2XYZ'])                  # what the first version did: one entry lost

    def test_domains_can_be_named_by_hand(self):
        self.assertEqual(rp.parse_domains('I=23-217, II=218-341'), [('I', 23, 217), ('II', 218, 341)])
        self.assertEqual(rp.parse_domains(''), [])
        with self.assertRaises(SystemExit):
            rp.parse_domains('I=23')


class Download(unittest.TestCase):
    def setUp(self):
        self._get = rp.get
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        rp.get = self._get

    def test_a_path_in_place_of_an_id_is_refused(self):
        rp.get = lambda *a, **k: self.fail('no request should be made')
        for bad in ('../x', '1ab', '12345', '', '1AB/', 'a123'):
            self.assertEqual(rp.download(bad, self.tmp)[1], False, bad)
        self.assertEqual(os.listdir(self.tmp), [])

    def test_an_error_page_is_not_a_structure(self):
        rp.get = lambda *a, **k: b'<html>' + b'x' * 5000 + b'</html>'
        self.assertEqual(rp.download('1ABC', self.tmp), ('1ABC', False, 'not an mmCIF file'))
        self.assertEqual(os.listdir(self.tmp), [])

    def test_a_good_file_is_written_whole_and_then_cached(self):
        rp.get = lambda *a, **k: b'data_1ABC\n' + b'x' * 3000
        self.assertEqual(rp.download('1ABC', self.tmp), ('1ABC', True, 'downloaded'))
        self.assertEqual(os.listdir(self.tmp), ['1ABC.cif'])             # no .part file left behind
        rp.get = lambda *a, **k: self.fail('second call must use the cache')
        self.assertEqual(rp.download('1ABC', self.tmp), ('1ABC', True, 'cached'))


if __name__ == '__main__':
    unittest.main()
