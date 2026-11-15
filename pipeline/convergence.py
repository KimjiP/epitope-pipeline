#!/usr/bin/env python3
"""Where do the binders of one antigen converge, and is that more than chance?

An epitope is made of residues that are far apart in sequence and close in space, so a stretch
of sequence (the first version used 24 residues) states it badly. Here a PATCH is every antigen
residue whose C-alpha lies within RADIUS of a centre residue.

A binder HITS a patch when at least MIN_RES of its footprint residues, and at least MIN_SHARE of
them, lie inside. The best centre is the one with the most hits. A footprint is the set of
residues the binder touches in at least half of its structures.

A count of hits means little alone: a footprint of 25 residues on a 600-residue protein lands
in some patch whatever it does. So the count is judged against a BASELINE. Every binder's
footprint is replaced by a compact footprint of the same size at a random place on the exposed
surface, and the best count is found again, many times. How often chance does as well is the p
value.

The random numbers come from a small generator with fixed arithmetic, so a replay gives the
same baseline on every Python version.

Stdlib only.
"""
import math

RADIUS = 10.0        # A, between C-alpha atoms: a sphere about 20 A across holds one epitope
MIN_RES = 3          # a binder needs at least this many footprint residues in the patch ...
MIN_SHARE = 0.25     # ... and at least this share of its footprint


class Rng:
    """64-bit linear congruential generator. Python's own `random` has changed between versions."""

    def __init__(self, seed):
        self.s = (seed * 6364136223846793005 + 1442695040888963407) & 0xFFFFFFFFFFFFFFFF

    def below(self, n):
        self.s = (self.s * 6364136223846793005 + 1442695040888963407) & 0xFFFFFFFFFFFFFFFF
        return (self.s >> 11) % n


def neighbours(ca, radius=RADIUS):
    """{residue: frozenset of residues whose C-alpha is within `radius` of its own}."""
    ids = sorted(ca)
    return {a: frozenset(b for b in ids if math.dist(ca[a], ca[b]) <= radius) for a in ids}


def hitters(sphere, footprints, min_res=MIN_RES, min_share=MIN_SHARE):
    """Names of the footprints that hit the sphere."""
    out = []
    for name, fp in footprints.items():
        k = len(fp & sphere)
        if k >= min_res and k >= min_share * len(fp):
            out.append(name)
    return out


def _usable(footprints, ca, min_res=MIN_RES):
    fps = {n: frozenset(r for r in fp if r in ca) for n, fp in footprints.items()}
    return {n: fp for n, fp in fps.items() if len(fp) >= min_res}


def find_patch(footprints, ca, nb=None):
    """The best patch for {name: set of residues}; None when no binder has enough residues.

    Ties go to the centre with the larger summed share, then to the lower residue number."""
    nb = nb or neighbours(ca)
    fps = _usable(footprints, ca)
    best = None
    for centre in sorted(nb):
        hit = hitters(nb[centre], fps)
        if not hit:
            continue
        key = (len(hit), sum(len(fps[n] & nb[centre]) / len(fps[n]) for n in hit))
        if best is None or key > best[0]:
            best = (key, centre, hit)
    if best is None:
        return None
    _, centre, hit = best
    touched = set()
    for n in hit:
        touched |= fps[n] & nb[centre]
    return {'centre': centre, 'hits': len(hit), 'considered': len(fps), 'binders': sorted(hit),
            'residues': sorted(nb[centre]), 'touched': sorted(touched)}


def baseline(footprints, ca, exposed, nb, draws=200, seed=7):
    """Best hit count, `draws` times, when every binder lands at a random place instead.

    Each random footprint is the k exposed residues nearest to a random exposed residue, k being
    the size of the real footprint: compact, the same size, and on the surface."""
    ids = sorted(r for r in exposed if r in ca)
    if not ids:
        return []
    order = {a: sorted(ids, key=lambda b: (math.dist(ca[a], ca[b]), b)) for a in ids}
    sizes = sorted(len(fp) for fp in _usable(footprints, ca).values())
    rng = Rng(seed)
    out = []
    for _ in range(draws):
        fake = {i: frozenset(order[ids[rng.below(len(ids))]][:k]) for i, k in enumerate(sizes)}
        out.append(max((len(hitters(nb[c], fake)) for c in nb), default=0))
    return out


def chance(null, observed):
    """Summary of the baseline against the observed hit count."""
    s = sorted(null)
    if not s:
        return None
    return {'draws': len(s), 'median': s[len(s) // 2], 'p95': s[int(0.95 * (len(s) - 1))],
            'p': round((1 + sum(1 for x in s if x >= observed)) / (1 + len(s)), 4)}


def pick_hotspots(candidates, counts, rsasa, density, min_binders=2, min_rsasa=0.20, n=8):
    """Residues to aim at: touched by at least `min_binders` binders and at least `min_rsasa`
    exposed in the free antigen, best first. Fewer than `n` is fine: a buried residue, or one no
    binder touches, is not a hotspot."""
    ok = [r for r in candidates if counts.get(r, 0) >= min_binders and rsasa.get(r, 0.0) >= min_rsasa]
    ok.sort(key=lambda r: (-counts[r], -density.get(r, 0), r))
    return sorted(ok[:n])


def runs_text(ids):
    """[5, 6, 7, 10] -> '5-7,10'"""
    ids = sorted(ids)
    if not ids:
        return ''
    out, st, pv = [], ids[0], ids[0]
    for r in ids[1:]:
        if r - pv > 1:
            out.append((st, pv)); st = r
        pv = r
    out.append((st, pv))
    return ','.join(f'{a}-{b}' if a != b else str(a) for a, b in out)
