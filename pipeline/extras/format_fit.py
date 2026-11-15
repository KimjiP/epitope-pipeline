#!/usr/bin/env python3
"""Measure the shape of an epitope or a candidate patch, for format choice.

Usage:
    python3 format_fit.py STRUCTURE CHAIN --set NAME=557,558,560 [--set ...]
                          [--offset N] [--json OUT.json]

Example:
    python3 format_fit.py 1N8Z.cif C --set trastuzumab=557,558,560,... --json ref.json

Why
---
A Fab and a VHH do not want the same surface. A Fab carries six loops across
two chains and can span a wide, broken-up epitope. A VHH carries three loops on
one domain, and its long CDR3 can enter a cleft a Fab cannot reach. So the
shape of a patch, not only its size, decides which format has a chance.

This measures the shape. It does not decide anything by itself: run it on
binders that already exist for your antigen, and compare.

What it computes per residue set
--------------------------------
    area         solvent-accessible surface of the set, antigen alone, A^2
    span         largest distance between any two residue centroids, A
    segments     contiguous stretches in sequence
    planarity    RMS distance of the centroids from their best-fit plane, A
                 (low = flat, high = corrugated or curved)
    convexity    centre minus rim, measured along the outward normal, A
                 (positive = the middle bulges out, negative = a dish or cleft)
    burial       mean count of antigen heavy atoms within 10 A of a centroid
                 (high = recessed, surrounded by its own protein)

`--offset N` is added to every residue number you pass, for structures that use
a different numbering convention from your notes.

Stdlib only.
"""
import sys, json, math, argparse
from collections import defaultdict

AA3 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
       'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
       'THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
RAD = {'C':1.70,'N':1.55,'O':1.52,'S':1.80,'SE':1.90}
HYDROPHOBIC = set('FWYLIVM')
PROBE, NPTS = 1.4, 96


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
    if not cols:
        sys.exit('no _atom_site loop found')
    ix = {c: n for n, c in enumerate(cols)}
    out = []
    for r in rows:
        if r[ix['group_PDB']] != 'ATOM' or r[ix['type_symbol']] == 'H':
            continue
        out.append(dict(ch=r[ix['auth_asym_id']], rs=int(r[ix['auth_seq_id']]),
                        aa=r[ix['label_comp_id']], at=r[ix['label_atom_id']],
                        el=r[ix['type_symbol']],
                        xyz=(float(r[ix['Cartn_x']]), float(r[ix['Cartn_y']]),
                             float(r[ix['Cartn_z']]))))
    return out


def unit_sphere(n):
    pts, ga = [], math.pi * (3 - 5 ** 0.5)
    for k in range(n):
        y = 1 - 2 * k / (n - 1)
        r = max(0.0, 1 - y * y) ** 0.5
        th = ga * k
        pts.append((math.cos(th) * r, y, math.sin(th) * r))
    return pts
SPH = unit_sphere(NPTS)


def sasa(atoms):
    xyz = [a['xyz'] for a in atoms]
    rad = [RAD.get(a['el'], 1.7) + PROBE for a in atoms]
    cell = 2 * max(rad) + 0.1
    grid = defaultdict(list)
    for n, (x, y, z) in enumerate(xyz):
        grid[(int(x//cell), int(y//cell), int(z//cell))].append(n)
    out = []
    for n, (p, rn) in enumerate(zip(xyz, rad)):
        cx, cy, cz = (int(v//cell) for v in p); nb = []
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for m in grid.get((cx+dx, cy+dy, cz+dz), ()):
                        if m != n and math.dist(p, xyz[m]) < rn + rad[m]:
                            nb.append(m)
        acc = 0
        for sx, sy, sz in SPH:
            q = (p[0] + sx*rn, p[1] + sy*rn, p[2] + sz*rn)
            if all(math.dist(q, xyz[m]) >= rad[m] for m in nb):
                acc += 1
        out.append(4 * math.pi * rn * rn * acc / NPTS)
    return out


def jacobi3(a):
    """Eigenvalues and eigenvectors of a symmetric 3x3 matrix."""
    a = [row[:] for row in a]
    v = [[1.0 if i == j else 0.0 for j in range(3)] for i in range(3)]
    for _ in range(64):
        off = sum(a[i][j] ** 2 for i in range(3) for j in range(3) if i != j)
        if off < 1e-18:
            break
        p, q = max(((i, j) for i in range(3) for j in range(3) if i < j),
                   key=lambda t: abs(a[t[0]][t[1]]))
        if abs(a[p][q]) < 1e-18:
            break
        theta = (a[q][q] - a[p][p]) / (2 * a[p][q])
        t = (1 if theta >= 0 else -1) / (abs(theta) + math.sqrt(theta * theta + 1))
        c = 1 / math.sqrt(t * t + 1); s = t * c
        for k in range(3):
            akp, akq = a[k][p], a[k][q]
            a[k][p] = c * akp - s * akq
            a[k][q] = s * akp + c * akq
        for k in range(3):
            apk, aqk = a[p][k], a[q][k]
            a[p][k] = c * apk - s * aqk
            a[q][k] = s * apk + c * aqk
        for k in range(3):
            vkp, vkq = v[k][p], v[k][q]
            v[k][p] = c * vkp - s * vkq
            v[k][q] = s * vkp + c * vkq
    vals = [a[i][i] for i in range(3)]
    vecs = [[v[0][i], v[1][i], v[2][i]] for i in range(3)]
    order = sorted(range(3), key=lambda i: vals[i])
    return [vals[i] for i in order], [vecs[i] for i in order]


def shape(cen_list, centre_all, cen_by_res):
    """Plane fit over residue centroids; returns planarity and convexity."""
    n = len(cen_list)
    px = sum(p[0] for p in cen_list) / n
    py = sum(p[1] for p in cen_list) / n
    pz = sum(p[2] for p in cen_list) / n
    P = (px, py, pz)
    cov = [[0.0] * 3 for _ in range(3)]
    for p in cen_list:
        d = (p[0] - px, p[1] - py, p[2] - pz)
        for i in range(3):
            for j in range(3):
                cov[i][j] += d[i] * d[j]
    for i in range(3):
        for j in range(3):
            cov[i][j] /= n
    vals, vecs = jacobi3(cov)
    normal = vecs[0]                       # smallest variance = plane normal
    out = (P[0] - centre_all[0], P[1] - centre_all[1], P[2] - centre_all[2])
    if sum(normal[i] * out[i] for i in range(3)) < 0:
        normal = [-normal[i] for i in range(3)]

    signed, lateral = [], []
    for p in cen_list:
        d = (p[0] - px, p[1] - py, p[2] - pz)
        s = sum(d[i] * normal[i] for i in range(3))
        signed.append(s)
        lat = math.sqrt(max(0.0, sum(x * x for x in d) - s * s))
        lateral.append(lat)
    planarity = math.sqrt(sum(s * s for s in signed) / n)
    rmax = max(lateral) or 1.0
    core = [s for s, l in zip(signed, lateral) if l <= 0.45 * rmax]
    rim = [s for s, l in zip(signed, lateral) if l >= 0.70 * rmax]
    convexity = ((sum(core) / len(core)) - (sum(rim) / len(rim))) if core and rim else 0.0
    return planarity, convexity, P


def runs(ids):
    out, st, pv = [], ids[0], ids[0]
    for r in ids[1:]:
        if r - pv > 1:
            out.append((st, pv)); st = r
        pv = r
    out.append((st, pv))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('structure')
    ap.add_argument('chain')
    ap.add_argument('--set', action='append', required=True, metavar='NAME=RESIDUES')
    ap.add_argument('--offset', type=int, default=0,
                    help='added to every residue number you pass')
    ap.add_argument('--json')
    args = ap.parse_args()

    atoms = parse(args.structure)
    ch = [a for a in atoms if a['ch'] == args.chain]
    if not ch:
        sys.exit(f'chain {args.chain} absent. present: {sorted({a["ch"] for a in atoms})}')
    for a, v in zip(ch, sasa(ch)):
        a['sasa'] = v

    res = defaultdict(list)
    for a in ch:
        res[a['rs']].append(a)
    aa1 = {r: AA3.get(res[r][0]['aa'], 'X') for r in res}

    def centroid(r):
        sc = [a for a in res[r] if a['at'] not in ('N', 'C', 'O')] or res[r]
        return tuple(sum(a['xyz'][k] for a in sc) / len(sc) for k in range(3))
    cen = {r: centroid(r) for r in res}

    allxyz = [a['xyz'] for a in ch]
    g = (sum(p[0] for p in allxyz) / len(allxyz),
         sum(p[1] for p in allxyz) / len(allxyz),
         sum(p[2] for p in allxyz) / len(allxyz))

    cell = 10.0
    grid = defaultdict(list)
    for n, (x, y, z) in enumerate(allxyz):
        grid[(int(x//cell), int(y//cell), int(z//cell))].append(n)

    def burial(p):
        cx, cy, cz = (int(v//cell) for v in p); c = 0
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for n in grid.get((cx+dx, cy+dy, cz+dz), ()):
                        if math.dist(p, allxyz[n]) <= 10.0:
                            c += 1
        return c

    results = []
    for spec in args.set:
        name, _, rl = spec.partition('=')
        ids = sorted({int(v) + args.offset for v in rl.replace(' ', '').split(',') if v})
        have = [r for r in ids if r in res]
        missing = [r for r in ids if r not in res]
        if len(have) < 4:
            print(f'! {name}: only {len(have)} of {len(ids)} residues are modelled; skipped')
            continue
        pts = [cen[r] for r in have]
        area = sum(sum(a['sasa'] for a in res[r]) for r in have)
        span = max(math.dist(p, q) for p in pts for q in pts)
        planarity, convexity, P = shape(pts, g, cen)
        bur = sum(burial(p) for p in pts) / len(pts)
        hyd = [f'{aa1[r]}{r}' for r in have if aa1[r] in HYDROPHOBIC]
        segs = runs(have)
        results.append({
            'name': name, 'structure': args.structure, 'chain': args.chain,
            'residues': have, 'missing': missing,
            'n_residues': len(have),
            'area': round(area, 1), 'span': round(span, 1),
            'segments': len(segs),
            'segment_ranges': [f'{a}-{b}' if a != b else str(a) for a, b in segs],
            'planarity': round(planarity, 2),
            'convexity': round(convexity, 2),
            'burial': round(bur, 1),
            'hydrophobic': hyd,
        })

    hdr = (f"{'set':<22} {'res':>4} {'area':>7} {'span':>6} {'segs':>5} "
           f"{'planar':>7} {'convex':>7} {'burial':>7}  hydrophobic")
    print(f"{args.structure} chain {args.chain}"
          + (f"  (offset {args.offset:+d})" if args.offset else ""))
    print(hdr); print('-' * len(hdr))
    for r in results:
        print(f"{r['name']:<22} {r['n_residues']:>4} {r['area']:>7.0f} {r['span']:>6.1f} "
              f"{r['segments']:>5} {r['planarity']:>7.2f} {r['convexity']:>7.2f} "
              f"{r['burial']:>7.1f}  {','.join(r['hydrophobic']) or '-'}")

    if args.json:
        with open(args.json, 'w') as fh:
            json.dump({'structure': args.structure, 'chain': args.chain,
                       'offset': args.offset, 'sets': results}, fh, indent=2)
            fh.write('\n')
        print(f"\nwrote {args.json}")


if __name__ == '__main__':
    main()
