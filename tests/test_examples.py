"""Scientific regression checks on the two saved examples (HER2, EGFR).

These read the files that ./reproduce.sh writes (.replay/<name>/), so they run after it and skip when
those files are absent. They check facts a structural biologist would check, not exact numbers (the
numbers are compared by reproduce.sh itself):

  * the approved drugs sit in the domains their papers report (Cho 2003; Franklin 2004; Li 2005)
  * one antibody is one binder, whatever it is called and however the entry was solved
  * every binder is drawn on the antigen, and every hotspot is exposed and touched

Run:  ./reproduce.sh && python3 -I -m unittest discover -s tests -v
"""
import json
import os
import unittest

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
WORK = os.environ.get('EPITOPE_WORK') or os.path.join(ROOT, '.replay')

# Domain borders in the mature numbering that the reference structure uses (UniProt minus the signal peptide)
HER2 = {'I': (1, 195), 'II': (196, 319), 'III': (320, 488), 'IV': (489, 607)}
EGFR = {'I': (1, 165), 'II': (166, 309), 'III': (310, 481), 'IV': (482, 621)}


def load(name):
    d = os.path.join(WORK, name)
    if not os.path.exists(os.path.join(d, 'atlas-out', 'atlas.json')):
        return None
    A = json.load(open(os.path.join(d, 'atlas-out', 'atlas.json')))
    S = json.load(open(os.path.join(d, 'summary.json')))
    P = json.load(open(os.path.join(d, 'props.json')))['residues']
    # Output left by an older version of the code is not what these checks are about.
    committed = os.path.join(ROOT, 'examples', name, 'summary.json')
    if A.get('schema', 0) < 3 or (os.path.exists(committed) and any(
            S.get(k) != json.load(open(committed)).get(k) for k in ('structures_used', 'units', 'distinct_binders'))):
        return 'stale'
    return A, S, P


def cluster_of(A, pdb, cdr3=None):
    """The cluster that holds a unit of `pdb` (with this heavy-chain CDR3, if given)."""
    for u in A['units']:
        if u['pdb'] != pdb:
            continue
        if cdr3 and not any(x[1] == cdr3 for x in (u.get('v') or {}).get('H', [])):
            continue
        return next(c for c in A['clusters'] if c['id'] == u['cluster'])
    raise AssertionError(f'no unit of {pdb}')


def share_in(residues, lo, hi):
    return sum(1 for r in residues if lo <= r <= hi) / max(1, len(residues))


class Common:
    NAME = None
    DOMAINS = None

    @classmethod
    def setUpClass(cls):
        got = load(cls.NAME)
        if got is None:
            raise unittest.SkipTest(f'no replay output for {cls.NAME}; run ./reproduce.sh first')
        if got == 'stale':
            raise unittest.SkipTest(f'the replay output for {cls.NAME} is from an older version; run ./reproduce.sh again')
        cls.A, cls.S, cls.P = got

    def test_every_binder_is_drawn_on_the_antigen(self):
        worst = max(u['rmsd'] for u in self.A['units'])
        self.assertLessEqual(worst, 3.0, 'a binder is placed with a fit error above 3 A')

    def test_every_hotspot_is_exposed_and_touched(self):
        off = self.S['reference_offset']
        counts = {int(k): v for k, v in self.A['counts_distinct'].items()}
        self.assertTrue(self.S['hotspots_uniprot'], 'no hotspots')
        for h in self.S['hotspots_uniprot']:
            r = int(h[1:]) - off
            self.assertGreaterEqual(self.P[str(r)]['rsasa'], 0.20, h)
            self.assertGreaterEqual(counts.get(r, 0), 2, h)

    def test_the_patch_comes_with_its_baseline(self):
        p = self.S['most_targeted_patch']
        self.assertIsNotNone(p['chance'])
        self.assertIn('median', p['chance'])
        self.assertGreaterEqual(p['of'], p['hits'])

    def test_counts_agree(self):
        self.assertEqual(self.S['distinct_binders'], len(self.A['clusters']))
        self.assertEqual(self.S['units'], len(self.A['units']))
        kinds = self.S['distinct_antibodies'] + self.S['distinct_scaffolds'] + self.S['distinct_other']
        self.assertEqual(kinds, self.S['distinct_binders'])

    def test_no_binder_is_a_relative_of_the_antigen(self):
        for u in self.A['units']:
            self.assertFalse(any(w in ' '.join(u['descriptions']).lower() for w in self.RELATIVES), u['name'])


class Her2(Common, unittest.TestCase):
    NAME = 'her2'
    RELATIVES = ('receptor tyrosine-protein kinase erbb-3', 'receptor tyrosine-protein kinase erbb-4',
                 'epidermal growth factor receptor')

    def test_trastuzumab_binds_domain_iv(self):
        c = cluster_of(self.A, '1N8Z')
        self.assertEqual(share_in(c['epitope_core'], *self.DOM('IV')), 1.0)

    def test_pertuzumab_binds_domain_ii(self):
        c = cluster_of(self.A, '1S78')
        self.assertGreaterEqual(share_in(c['epitope_core'], *self.DOM('II')), 0.9)

    def test_low_resolution_entries_are_left_out(self):
        # 8FFJ (zanidatamab, 7.5 A), 6OGE (4.36 A) and 5KWG (4.3 A) are skipped by the 4 A cut.
        # The rule that merges zanidatamab's second arm into trastuzumab is tested on the real
        # sequences in test_binders.py.
        self.assertGreaterEqual(self.S['entries_skipped'].get('resolution worse than 4 A', 0), 3)
        self.assertLessEqual(max(u['resolution'] or 0 for u in self.A['units']), 4.0)
        self.assertNotIn('8FFJ', {u['pdb'] for u in self.A['units']})

    def test_bh1_is_not_trastuzumab(self):
        # 3BE1: trastuzumab's heavy chain with a light chain 7% different
        self.assertNotEqual(cluster_of(self.A, '3BE1')['id'], cluster_of(self.A, '1N8Z')['id'])

    def test_trastuzumab_is_one_binder_over_all_its_structures(self):
        self.assertGreaterEqual(len(cluster_of(self.A, '1N8Z')['pdbs']), 5)

    @staticmethod
    def DOM(name):
        return HER2[name]


class Egfr(Common, unittest.TestCase):
    NAME = 'egfr'
    RELATIVES = ('receptor tyrosine-protein kinase erbb-2', 'receptor tyrosine-protein kinase erbb-3')

    def test_cetuximab_binds_domain_iii_only(self):
        c = cluster_of(self.A, '1YY9')
        self.assertEqual(share_in(c['epitope_core'], *EGFR['III']), 1.0)

    def test_imc_11f8_and_necitumumab_are_one_antibody(self):
        a, b, c = (cluster_of(self.A, p)['id'] for p in ('3B2U', '3B2V', '6B3S'))
        self.assertEqual(a, b)
        self.assertEqual(b, c)

    def test_the_natural_growth_factors_are_not_counted_as_antibodies(self):
        egf = cluster_of(self.A, '1IVO')
        self.assertEqual(egf['kind'], 'other')
        self.assertNotEqual(egf['id'], cluster_of(self.A, '1YY9')['id'])

    def test_antibodies_converge_on_domain_iii(self):
        p = self.S['most_targeted_patch']
        self.assertEqual(p['basis'], 'antibodies and nanobodies')
        self.assertEqual(p['domain'], 'III')
        self.assertLess(p['chance']['p'], 0.05)

    def test_the_example_reads_every_extracellular_structure(self):
        self.assertFalse(self.S['sampled'])
        self.assertEqual(self.S['structures_in_region'], self.S['structures_used'] + len(self.S['quality']['downloads_failed']))


if __name__ == '__main__':
    unittest.main()
