"""Tests for the antibody logic in atlas.py: finding variable domains, comparing them, and deciding
which binding units are one binder. The sequences are real (see sequences.py).

Run:  python3 -I -m unittest discover -s tests -v
"""
import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'pipeline'))
sys.path.insert(0, HERE)

import atlas                       # noqa: E402
import sequences as S              # noqa: E402


def mutate(seq, positions, to='A'):
    s = list(seq)
    for p in positions:
        s[p] = to if s[p] != to else 'G'
    return ''.join(s)


def unit(uid, h=None, l=None, fmt='Fab', kind='antibody', seqs=None):
    """A binding unit as atlas.py builds it, with only the fields the comparison reads."""
    v = {}
    if h:
        v['H'] = [h]
    if l:
        v['L'] = [l]
    return {'id': uid, 'kind': kind, 'format': fmt, 'v': v or None, 'seqs': seqs or []}


TZ = (S.TZ_H, S.TZ_H_CDR3), (S.TZ_L, S.TZ_L_CDR3)
BH1 = (S.BH1_H, S.BH1_H_CDR3), (S.BH1_L, S.BH1_L_CDR3)
PZ = (S.PZ_H, S.PZ_H_CDR3), (S.PZ_L, S.PZ_L_CDR3)
IMC = (S.IMC_H, S.IMC_H_CDR3), (S.IMC_L, S.IMC_L_CDR3)
NEC = (S.NEC_H, S.NEC_H_CDR3), (S.NEC_L, S.NEC_L_CDR3)


class FindDomains(unittest.TestCase):
    def test_heavy_and_light_found_in_whole_chains(self):
        for vseq, cdr3, kind, const in ((S.TZ_H, S.TZ_H_CDR3, 'H', S.CH1),
                                        (S.TZ_L, S.TZ_L_CDR3, 'L', S.CKAPPA),
                                        (S.PZ_H, S.PZ_H_CDR3, 'H', S.CH1),
                                        (S.IMC_L, S.IMC_L_CDR3, 'L', S.CKAPPA)):
            d = atlas.v_domains(vseq + const)
            self.assertEqual(len(d), 1)
            self.assertEqual((d[0][0], d[0][1], d[0][2]), (kind, vseq, cdr3))

    def test_nanobody_is_a_heavy_type_domain(self):
        d = atlas.v_domains(S.VHH7D12)
        self.assertEqual([x[0] for x in d], ['H'])

    def test_the_last_residues_of_the_j_region_may_be_missing(self):
        # three Fab heavy chains of the HER2 set end at ...LVTVS or ...LVTVF in the model
        d = atlas.v_domains(S.TZ_H[:-1])
        self.assertEqual(len(d), 1)
        d = atlas.v_domains(S.TZ_H[:-2] + 'F')
        self.assertEqual(len(d), 1)

    def test_an_unusual_light_chain_j_region(self):
        # 3WSQ: ...FGSGTFLTVTVAAPS, where a lambda chain has FGGGTKLTVL
        vl = S.TZ_L[:-10] + 'FGSGTFLTVT'
        self.assertEqual([x[0] for x in atlas.v_domains(vl + 'VAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKV')], ['L'])

    def test_light_chains_with_odd_j_regions(self):
        # 6J14 (mouse, ...FGGGTKLEMK) and cemiplimab (7WVM, ...FGPGTVVDFR): both missed at first
        for j in ('FGGGTKLEMK', 'FGPGTVVDFR', 'FGGGTKLTVL', 'FGQGTRLEIK', 'FGAGTKLELK'):
            vl = S.TZ_L[:-10] + j
            self.assertEqual([x[0] for x in atlas.v_domains(vl + S.CKAPPA)], ['L'], j)

    def test_an_scfv_has_two_domains_and_the_linker_is_skipped(self):
        d = atlas.v_domains(S.TZ_H + 'GGGGSGGGGSGGGGS' + S.TZ_L)
        self.assertEqual([x[0] for x in d], ['H', 'L'])
        self.assertEqual(d[1][1], S.TZ_L)

    def test_constant_domains_and_other_proteins_have_none(self):
        self.assertEqual(atlas.v_domains(S.CH1 + S.CKAPPA), [])
        self.assertEqual(atlas.v_domains(S.EGF), [])
        self.assertEqual(atlas.v_domains(''), [])


class Identity(unittest.TestCase):
    def test_overlap_identity(self):
        a = S.TZ_H
        self.assertEqual(atlas.overlap_identity(a, a), 1.0)
        self.assertEqual(atlas.overlap_identity(a, a[4:]), 1.0)               # starts later
        self.assertEqual(atlas.overlap_identity(a[:-3], a), 1.0)              # stops sooner
        self.assertAlmostEqual(atlas.overlap_identity(a, mutate(a, [10])), (len(a) - 1) / len(a))
        self.assertLess(atlas.overlap_identity(S.TZ_H, S.PZ_H), 0.8)          # two different clones
        self.assertEqual(atlas.overlap_identity('', a), 0.0)

    def test_nw_identity_still_punishes_a_truncated_copy(self):
        a = S.TZ_H
        self.assertAlmostEqual(atlas.nw_identity(a, a[:-12]), (len(a) - 12) / len(a))


class Classify(unittest.TestCase):
    def test_kind_and_format(self):
        c = atlas.classify
        both = {'H': [('x', '')], 'L': [('y', '')]}
        self.assertEqual(c(['Fab heavy chain', 'Fab light chain'], '', [220, 214], both), ('antibody', 'Fab'))
        self.assertEqual(c(['heavy', 'light'], '', [118, 108], both), ('antibody', 'Fv'))
        self.assertEqual(c(['Antibody'], '', [220, 214], both, one_chain_both=True), ('antibody', 'scFv'))
        self.assertEqual(c(['Nanobody'], '', [124], {'H': [('x', '')]}), ('antibody', 'VHH'))
        self.assertEqual(c(['heavy chain'], '', [221], {'H': [('x', '')]}), ('antibody', 'antibody fragment'))
        # no variable domain found (the J region is not modelled): the name decides
        self.assertEqual(c(['anti-TNF nanobodies, derived from Ozoralizumab'], '', [224]), ('antibody', 'VHH'))

    def test_the_sequence_beats_the_name(self):
        # the PDB name says only 'Antibody' or 'heavy chain'; the variable domain settles it
        self.assertEqual(atlas.classify(['light chain of something'], '', [214], {'L': [('y', '')]})[0], 'antibody')
        self.assertEqual(atlas.classify(['Protein X'], '', [214], {'L': [('y', '')]})[0], 'antibody')

    def test_scaffolds_and_other_proteins(self):
        c = atlas.classify
        self.assertEqual(c(['Designed Ankyrin Repeat Protein 9_29'], '', [160]), ('scaffold', 'DARPin'))
        self.assertEqual(c(['Immunoglobulin G-binding protein A'], '', [58]), ('scaffold', 'Affibody'))
        self.assertEqual(c(['Ig gamma-1 chain C region'], '', [220]), ('scaffold', 'Fcab (engineered Fc)'))
        # the natural ligand of EGFR is not a scaffold, and no longer labelled one
        self.assertEqual(c(['Epidermal growth factor'], '', [53]), ('other', 'other protein'))


class Clusters(unittest.TestCase):
    def groups(self, units, thr=0.95):
        g, _ = atlas.cluster_units(units, thr)
        return sorted(sorted(x) for x in g)

    def test_one_drug_in_three_formats_is_one_binder(self):
        us = [unit('U01', *TZ),
              unit('U02', (mutate(S.TZ_H, [5, 50]), S.TZ_H_CDR3), (mutate(S.TZ_L, [7]), S.TZ_L_CDR3), fmt='Fv'),
              unit('U03', (S.TZ_H, S.TZ_H_CDR3), fmt='antibody fragment')]          # a VH with no modelled VL
        self.assertEqual(self.groups(us), [['U01', 'U02', 'U03']])

    def test_the_same_drug_under_two_names_is_one_binder(self):
        # IMC-11F8 (3B2U, 3B2V) and necitumumab (6B3S): the first version counted three binders
        self.assertEqual(self.groups([unit('U01', *IMC), unit('U02', *NEC)]), [['U01', 'U02']])

    def test_a_different_clone_is_a_different_binder(self):
        us = [unit('U01', *TZ), unit('U02', *PZ), unit('U03', *IMC)]
        self.assertEqual(self.groups(us), [['U01'], ['U02'], ['U03']])

    def test_a_shared_heavy_chain_does_not_hide_a_different_light_chain(self):
        # bH1: trastuzumab's heavy chain, a light chain 7% different. The mean of the two
        # chains would call it trastuzumab; the lower of the two does not.
        self.assertEqual(self.groups([unit('U01', *TZ), unit('U02', *BH1)]), [['U01'], ['U02']])

    def test_a_partial_unit_never_joins_two_clones_together(self):
        # a VH identical to the heavy chain of BOTH trastuzumab and bH1
        us = [unit('U01', *TZ), unit('U02', *BH1), unit('U03', (S.TZ_H, S.TZ_H_CDR3))]
        g = self.groups(us)
        self.assertEqual(len(g), 2)
        self.assertEqual(sum(len(x) for x in g), 3)
        self.assertNotEqual({u for x in g for u in x if 'U01' in x}, {u for x in g for u in x if 'U02' in x})
        self.assertTrue(any(x == ['U01', 'U03'] for x in g))      # it goes to the larger group, by the rule

    def test_a_partial_unit_needs_the_same_cdr3(self):
        us = [unit('U01', *TZ), unit('U02', (S.TZ_H, 'SRAAAAAAAAAAY'))]
        self.assertEqual(self.groups(us), [['U01'], ['U02']])

    def test_nanobodies_group_among_themselves(self):
        v = (S.VHH7D12, S.VHH7D12_CDR3)
        us = [unit('U01', *TZ), unit('U02', v, fmt='VHH'), unit('U03', (mutate(S.VHH7D12, [20]), S.VHH7D12_CDR3), fmt='VHH')]
        self.assertEqual(self.groups(us), [['U01'], ['U02', 'U03']])

    def test_other_proteins_are_the_same_binder_whatever_the_construct(self):
        egf = unit('U01', kind='other', fmt='other protein', seqs=[S.EGF])
        short = unit('U02', kind='other', fmt='other protein', seqs=[S.EGF[2:-2]])       # modelled a little shorter
        other = unit('U03', kind='other', fmt='other protein', seqs=[S.CKAPPA[:60]])
        self.assertEqual(self.groups([egf, short, other]), [['U01', 'U02'], ['U03']])

    def test_an_antibody_never_matches_a_non_antibody(self):
        self.assertEqual(atlas.unit_identity(unit('U01', *TZ), unit('U02', kind='other', fmt='other protein', seqs=[S.EGF])), 0.0)

    def test_two_different_fcabs_stay_apart(self):
        # the Fc is the same everywhere; only the last 110 residues (CH3) carry the engineered loops
        base = (S.CH1 + S.CKAPPA)[:130]
        a = base + S.TZ_H
        b = base + S.PZ_H[:len(S.TZ_H)]
        fa = unit('U01', kind='scaffold', fmt='Fcab (engineered Fc)', seqs=[a])
        fb = unit('U02', kind='scaffold', fmt='Fcab (engineered Fc)', seqs=[b])
        self.assertLess(atlas.unit_identity(fa, fb), 0.95)
        self.assertGreaterEqual(atlas.unit_identity(fa, dict(fa, id='U03')), 0.95)


class Fit(unittest.TestCase):
    def test_a_domain_that_moved_is_fitted_locally(self):
        # a rigid body of 80 points plus a second body that sits 30 A away in the other structure
        import math
        body = [(float(i % 5) * 3.8, float((i // 5) % 4) * 3.8, float(i // 20) * 3.8) for i in range(80)]
        other = [(x + 60.0, y, z) for x, y, z in body]
        ref = {i + 1: p for i, p in enumerate(body + other)}
        moved = {i + 1: p for i, p in enumerate(body + [(x + 30.0, y + 10.0, z) for x, y, z in other])}
        epitope = [1, 2, 3, 6, 7]                       # on the first body
        R, t, rmsd, n, scope = atlas.local_fit(moved, ref, epitope, 15.0, ((), (), 99.0))
        self.assertEqual(scope, 'local')
        self.assertLess(rmsd, 0.01)
        self.assertGreater(n, 12)
        # the whole-chain fit is bad
        pts_m, pts_r = list(moved.values()), list(ref.values())
        self.assertGreater(atlas.superpose(pts_m, pts_r)[2], 5.0)

    def test_too_few_residues_falls_back_to_the_whole_chain_fit(self):
        ref = {i: (float(i), 0.0, 0.0) for i in range(1, 6)}
        res = atlas.local_fit(ref, ref, [1], 15.0, ('R', 't', 1.5))
        self.assertEqual(res, ('R', 't', 1.5, 5, 'global'))


if __name__ == '__main__':
    unittest.main()
