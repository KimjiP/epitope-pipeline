#!/usr/bin/env python3
"""Pick epitope sites by combining every evidence stream, then name anchors.

Usage:
    python3 pick_sites.py STRUCTURE CHAIN REF_CHAINS --props P.json --atlas A.json
                          [--window LO-HI] [--n 2] [--json OUT.json]

Example:
    python3 pick_sites.py 1N8Z.cif C AB --props ../viewer-atlas/props.json \
        --atlas ../viewer-atlas/atlas.json --window 490-607 --n 2

What goes into the score
------------------------
Per residue, five streams, each one measured elsewhere in this project:

    competition  in the reference antibody's footprint, or near it in space.
                 The challenge requires the design to block trastuzumab, so
                 this is a requirement and not a preference.
    precedent    how many distinct published binders touch it (counts_distinct in the atlas).
    exposure     relative SASA. Buried surface cannot be bound.
    anchor       exposed hydrophobic area, rSASA x Kyte-Doolittle where
                 positive. This is NOT evidence that antibodies go there —
                 the correlation with binder density is small (see surface_props.py).
                 It is included because a designed binder still needs
                 something to bury, which is a different claim.
    penalties    glycan sequon nearby, cysteine, unmodelled neighbours.

The weights are printed with the result and are arguments, not secrets. Change
them and the answer changes; that is the point of showing them.

Stdlib only.
"""
import sys, json, math, argparse
from collections import defaultdict

AA3 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
       'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
       'THR':'T','TRP':'W','TYR':'Y','VAL':'V','MSE':'M'}
HYDRO = set('FWYLIVM')

W = {'competition': 3.0, 'near_footprint': 1.2, 'precedent': 0.45,
     'exposure': 2.0, 'anchor': 1.6, 'burial': -1.2,
     'glycan': -3.0, 'cysteine': -1.0, 'edge': -2.5}


def parse(path):
    lines = open(path).read().splitlines()
    cols, rows, i = [], [], 0
    while i < len(lines):
        if lines[i].startswith('loop_'):
            j, tags = i + 1, []
            while j < len(lines) and lines[j].startswith('_'):
                tags.append(lines[j].strip()); j += 1
            if tags and tags[0].startswith('_atom_site.'):
                cols = [t.split('.', 1)[1] for t in tags]; k = j
                while k < len(lines) and lines[k].strip() and not lines[k].startswith('#'):
                    p = lines[k].split()
                    if len(p) == len(cols): rows.append(p)
                    k += 1
                i = k; continue
            i = j; continue
        i += 1
    ix = {c: n for n, c in enumerate(cols)}
    out = []
    for r in rows:
        if r[ix['group_PDB']] != 'ATOM' or r[ix['type_symbol']] == 'H':
            continue
        out.append(dict(ch=r[ix['auth_asym_id']], rs=int(r[ix['auth_seq_id']]),
                        aa=r[ix['label_comp_id']], at=r[ix['label_atom_id']],
                        xyz=(float(r[ix['Cartn_x']]), float(r[ix['Cartn_y']]),
                             float(r[ix['Cartn_z']]))))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('structure'); ap.add_argument('chain'); ap.add_argument('ref')
    ap.add_argument('--props', required=True); ap.add_argument('--atlas', required=True)
    ap.add_argument('--window', default='490-607')
    ap.add_argument('--n', type=int, default=2)
    ap.add_argument('--radius', type=float, default=11.0)
    ap.add_argument('--json')
    args = ap.parse_args()
    LO, HI = (int(v) for v in args.window.split('-'))

    P = json.load(open(args.props))['residues']
    A = json.load(open(args.atlas))
    cd = {int(k): v for k, v in A['counts_distinct'].items()}
    foot = (A.get('trastuzumab_footprint')
            or (A.get('footprints', {}).get('trastuzumab') or {}).get('union'))
    if not foot:
        sys.exit('this atlas holds no trastuzumab footprint. Build one with atlas.py ... --name trastuzumab=1N8Z')
    FOOT = set(foot)
    NDIST = A.get('n_distinct', 0)

    atoms = parse(args.structure)
    ag = [a for a in atoms if a['ch'] == args.chain]
    ref = [a for a in atoms if a['ch'] in list(args.ref)]
    res = defaultdict(list)
    for a in ag:
        res[a['rs']].append(a)
    ids = sorted(res)
    aa1 = {r: AA3.get(res[r][0]['aa'], 'X') for r in ids}

    def centroid(r):
        sc = [a for a in res[r] if a['at'] not in ('N', 'C', 'O')] or res[r]
        return tuple(sum(a['xyz'][k] for a in sc) / len(sc) for k in range(3))
    cen = {r: centroid(r) for r in ids}

    # atom pairs with the reference antibody, the energetic proxy
    pairs = defaultdict(int)
    if ref:
        cell = 4.5
        grid = defaultdict(list)
        for n, a in enumerate(ref):
            x, y, z = a['xyz']
            grid[(int(x//cell), int(y//cell), int(z//cell))].append(n)
        for a in ag:
            x, y, z = a['xyz']
            cx, cy, cz = int(x//cell), int(y//cell), int(z//cell)
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    for dz in (-1, 0, 1):
                        for n in grid.get((cx+dx, cy+dy, cz+dz), ()):
                            if math.dist(a['xyz'], ref[n]['xyz']) <= 4.5:
                                pairs[a['rs']] += 1

    # distance from each residue to the nearest footprint residue
    fcen = [cen[r] for r in FOOT if r in cen]
    dfoot = {r: (min(math.dist(cen[r], f) for f in fcen) if fcen else 99.0) for r in ids}

    # Residues beside a chain break or a terminus look far more exposed than
    # they are: the neighbouring protein simply is not in the model. 1N8Z chain
    # C stops at 607 although domain IV runs past it, and it is missing 581-590,
    # so 580, 591, 605, 606 and 607 all read as unusually open. Flag them.
    edges = set()
    for a, b in zip(ids, ids[1:]):
        if b - a > 1:
            edges.update([a, a - 1, b, b + 1])
    edges.update([ids[0], ids[0] + 1, ids[-1], ids[-1] - 1, ids[-1] - 2])
    edges = {r for r in edges if r in cen}
    print(f'exposure inflated by a chain break or terminus: '
          f'{", ".join(f"{aa1[r]}{r}" for r in sorted(edges) if LO <= r <= HI)}')

    glyc = [r for r in ids if aa1[r] == 'N' and aa1.get(r+1) != 'P' and aa1.get(r+2) in ('S', 'T')]
    gdist = {r: (min(math.dist(cen[r], cen[g]) for g in glyc) if glyc else 99.0) for r in ids}
    burials = [P[str(r)]['burial'] for r in ids if str(r) in P]
    bmin, bmax = min(burials), max(burials)

    def score(r):
        p = P.get(str(r))
        if not p:
            return None
        s = 0.0; why = []
        if r in FOOT:
            s += W['competition']; why.append('in footprint')
        elif dfoot[r] <= 9.0:
            s += W['near_footprint'] * (1 - dfoot[r] / 9.0); why.append(f'{dfoot[r]:.0f}A from it')
        s += W['precedent'] * min(cd.get(r, 0), 4)
        s += W['exposure'] * p['rsasa']
        s += W['anchor'] * p['hydrophobic']
        s += W['burial'] * ((p['burial'] - bmin) / (bmax - bmin))
        if gdist[r] <= 6.0:
            s += W['glycan']; why.append('glycan')
        if aa1[r] == 'C':
            s += W['cysteine']
        if r in edges:
            s += W['edge']; why.append('edge')
        return s

    scored = {r: score(r) for r in ids if LO <= r <= HI and score(r) is not None}
    exposed = {r: v for r, v in scored.items() if P[str(r)]['rsasa'] >= 0.20}

    print(f'{args.structure} chain {args.chain}, window {LO}-{HI}')
    print(f'reference antibody chains {args.ref}; {len(FOOT)} footprint residues; '
          f'{NDIST} distinct binders in the atlas')
    print('weights: ' + ', '.join(f'{k} {v:+g}' for k, v in W.items()))
    print()

    picked, used, seedbar = [], set(), set()
    for _ in range(args.n):
        cands = [r for r in exposed if r not in seedbar]
        if not cands:
            break
        best, bestpatch, bestval = None, None, -1e9
        for seed in cands:
            patch = [r for r in exposed
                     if r not in used and math.dist(cen[seed], cen[r]) <= args.radius]
            if len(patch) < 5:
                continue
            v = sum(exposed[r] for r in patch) / (len(patch) ** 0.35)
            if v > bestval:
                best, bestpatch, bestval = seed, patch, v
        if best is None:
            break
        picked.append((best, sorted(bestpatch), bestval))
        # Take out this patch, and bar any later seed from sitting on top of it,
        # but do not clear a halo: that starved the second site down to five
        # residues last time. Two epitopes may abut; they may not coincide.
        used.update(bestpatch)
        for r2 in list(exposed):
            if math.dist(cen[best], cen[r2]) <= 12.0:
                seedbar.add(r2)

    out = []
    for n, (seed, patch, v) in enumerate(picked, 1):
        inf = [r for r in patch if r in FOOT]
        anchors = sorted(patch, key=lambda r: -(
            2.0 * P[str(r)]['hydrophobic'] + 1.2 * P[str(r)]['rsasa']
            + 0.06 * pairs.get(r, 0) + (0.8 if r in FOOT else 0)))
        anchors = [r for r in anchors
                   if P[str(r)]['rsasa'] >= 0.30 and aa1[r] != 'C' and r not in edges][:5]
        runs, st, pv = [], patch[0], patch[0]
        for r in patch[1:]:
            if r - pv > 1:
                runs.append((st, pv)); st = r
            pv = r
        runs.append((st, pv))
        segs = ','.join(f'{a}-{b}' if a != b else str(a) for a, b in runs)
        print(f'=== SITE {n}  seed {aa1[seed]}{seed}  score {v:.1f} ===')
        print(f'  residues ({len(patch)}): {segs}')
        print(f'  in trastuzumab footprint: {len(inf)} of {len(patch)}  '
              f'[{", ".join(f"{aa1[r]}{r}" for r in inf)}]')
        print(f'  distinct binders touching any residue: '
              f'{max([cd.get(r, 0) for r in patch] + [0])} of {NDIST}')
        print(f'  mean rSASA {sum(P[str(r)]["rsasa"] for r in patch)/len(patch):.2f}  '
              f'mean burial {sum(P[str(r)]["burial"] for r in patch)/len(patch):.0f}  '
              f'exposed area {sum(P[str(r)]["area"] for r in patch):.0f} A2')
        print(f'  nearest glycan sequon: {min(gdist[r] for r in patch):.1f} A')
        print(f'  {"ANCHORS":>9}  {"rSASA":>6} {"area":>6} {"KD":>5} {"Fab pairs":>10} '
              f'{"binders":>8}  role')
        arec = []
        for r in anchors:
            p = P[str(r)]
            role = []
            if r in FOOT: role.append('trastuzumab contact')
            if aa1[r] in HYDRO and p['rsasa'] >= 0.35: role.append('exposed hydrophobic')
            if cd.get(r, 0) == 0: role.append('untouched by any binder')
            print(f'  {aa1[r] + str(r):>9}  {p["rsasa"]:>6.2f} {p["area"]:>6.0f} '
                  f'{p["kd"]:>5.1f} {pairs.get(r, 0):>10} {cd.get(r, 0):>8}  '
                  f'{"; ".join(role) or "support"}')
            arec.append({'residue': f'{aa1[r]}{r}', 'resi': r, 'aa': aa1[r],
                         'rsasa': p['rsasa'], 'area': p['area'], 'kd': p['kd'],
                         'fab_pairs': pairs.get(r, 0), 'binders': cd.get(r, 0),
                         'in_footprint': r in FOOT, 'role': role})
        print()
        out.append({'site': n, 'seed': f'{aa1[seed]}{seed}', 'score': round(v, 2),
                    'residues': patch, 'segments': segs,
                    'in_footprint': [f'{aa1[r]}{r}' for r in inf],
                    'anchors': arec,
                    'mean_rsasa': round(sum(P[str(r)]['rsasa'] for r in patch)/len(patch), 3),
                    'mean_burial': round(sum(P[str(r)]['burial'] for r in patch)/len(patch), 1),
                    'area': round(sum(P[str(r)]['area'] for r in patch), 1),
                    'nearest_glycan': round(min(gdist[r] for r in patch), 1)})

    if args.json:
        with open(args.json, 'w') as fh:
            json.dump({'structure': args.structure, 'window': [LO, HI],
                       'weights': W, 'radius': args.radius, 'sites': out}, fh, indent=2)
            fh.write('\n')
        print(f'wrote {args.json}')


if __name__ == '__main__':
    main()
