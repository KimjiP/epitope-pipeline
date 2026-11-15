"""Tests for convergence.py: the patch, the baseline, the hotspot filter.

Run:  python3 -I -m unittest discover -s tests -v
"""
import math
import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'pipeline'))

import convergence as cv           # noqa: E402


def grid(nx=10, ny=10, nz=6, step=5.0):
    """A block of 600 'residues' on a lattice, numbered from 1. At 5 A a 10 A sphere holds about
    30 of them, as in a real protein."""
    return {1 + i + nx * (j + ny * k): (i * step, j * step, k * step)
            for k in range(nz) for j in range(ny) for i in range(nx)}


def compact(ca, centre, size):
    """The `size` residues nearest to `centre`: a compact footprint."""
    return set(sorted(ca, key=lambda r: (math.dist(ca[centre], ca[r]), r))[:size])


class Patch(unittest.TestCase):
    def setUp(self):
        self.ca = grid()
        self.nb = cv.neighbours(self.ca)
        self.rng = cv.Rng(11)

    def footprints_around(self, centre, n, size_lo=12, size_hi=18):
        out = {}
        for i in range(n):
            size = size_lo + self.rng.below(size_hi - size_lo + 1)
            # each binder sits a little off the centre, as real ones do
            near = sorted(self.nb[centre])
            c = near[self.rng.below(len(near))]
            out[f'b{centre}_{i}'] = compact(self.ca, c, size)
        return out

    def test_a_planted_patch_is_found_and_beats_chance(self):
        centre = 1 + 5 + 10 * (5 + 10 * 3)             # near the middle of the block
        fp = self.footprints_around(centre, 9)
        for i in range(5):                              # five binders elsewhere
            fp[f'far{i}'] = compact(self.ca, self.rng.below(len(self.ca)) + 1, 14)
        res = cv.find_patch(fp, self.ca, self.nb)
        self.assertGreaterEqual(res['hits'], 9)
        self.assertLess(math.dist(self.ca[res['centre']], self.ca[centre]), 10.0)
        null = cv.baseline(fp, self.ca, set(self.ca), self.nb, draws=100)
        ch = cv.chance(null, res['hits'])
        self.assertLess(ch['p'], 0.05)
        self.assertLess(ch['median'], res['hits'])

    def test_binders_spread_over_the_antigen_are_not_above_chance(self):
        fp = {f'b{i}': compact(self.ca, self.rng.below(len(self.ca)) + 1, 15) for i in range(12)}
        res = cv.find_patch(fp, self.ca, self.nb)
        null = cv.baseline(fp, self.ca, set(self.ca), self.nb, draws=100)
        self.assertGreater(cv.chance(null, res['hits'])['p'], 0.05)

    def test_a_binder_that_grazes_the_edge_is_not_a_hit(self):
        centre = 1 + 5 + 10 * (5 + 10 * 3)
        sphere = self.nb[centre]
        edge = sorted(sphere, key=lambda r: -math.dist(self.ca[r], self.ca[centre]))[:2]   # 2 residues inside
        far = compact(self.ca, 1, 25)
        fp = {'graze': set(edge) | (far - sphere)}
        self.assertEqual(cv.hitters(sphere, fp), [])
        # three residues but under a quarter of a big footprint: still not a hit
        big = set(sorted(sphere)[:3]) | (compact(self.ca, 1, 60) - sphere)
        self.assertEqual(cv.hitters(sphere, {'big': big}), [])

    def test_too_small_a_footprint_is_ignored(self):
        self.assertIsNone(cv.find_patch({'tiny': {1, 2}}, self.ca, self.nb))


class Baseline(unittest.TestCase):
    def test_the_generator_gives_the_same_numbers_on_every_python(self):
        r = cv.Rng(7)
        self.assertEqual([r.below(1000) for _ in range(5)], [r2 for r2 in EXPECTED_RNG])

    def test_the_baseline_is_reproducible(self):
        ca = grid()
        nb = cv.neighbours(ca)
        fp = {f'b{i}': compact(ca, 5 + 9 * i, 14) for i in range(6)}
        a = cv.baseline(fp, ca, set(ca), nb, draws=20, seed=3)
        b = cv.baseline(fp, ca, set(ca), nb, draws=20, seed=3)
        self.assertEqual(a, b)

    def test_chance_summary(self):
        ch = cv.chance([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], 8)
        self.assertEqual(ch['draws'], 10)
        self.assertEqual(ch['p'], round((1 + 3) / 11, 4))
        self.assertIsNone(cv.chance([], 3))


EXPECTED_RNG = [213, 67, 337, 367, 433]          # computed once; the arithmetic is plain integer maths


class Hotspots(unittest.TestCase):
    def test_a_hotspot_is_touched_by_two_binders_and_exposed(self):
        counts = {10: 5, 11: 4, 12: 1, 13: 6, 14: 3}
        rsasa = {10: 0.5, 11: 0.05, 12: 0.9, 13: 0.25, 14: 0.2}
        got = cv.pick_hotspots([10, 11, 12, 13, 14, 15], counts, rsasa, {}, n=8)
        self.assertEqual(got, [10, 13, 14])         # 11 is buried, 12 has one binder, 15 has none

    def test_fewer_than_n_is_allowed_and_the_best_come_first(self):
        counts = {r: r for r in range(1, 20)}
        rsasa = {r: 0.6 for r in range(1, 20)}
        self.assertEqual(cv.pick_hotspots(range(1, 20), counts, rsasa, {}, n=3), [17, 18, 19])
        self.assertEqual(cv.pick_hotspots([5], {5: 1}, {5: 0.9}, {}), [])

    def test_runs_text(self):
        self.assertEqual(cv.runs_text([5, 6, 7, 10, 12, 13]), '5-7,10,12-13')
        self.assertEqual(cv.runs_text([]), '')


if __name__ == '__main__':
    unittest.main()
