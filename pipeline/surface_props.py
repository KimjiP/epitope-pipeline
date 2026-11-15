#!/usr/bin/env python3
"""Per-residue surface properties of an antigen, and whether they predict binders.

Usage:
    python3 surface_props.py STRUCTURE CHAIN [--atlas ATLAS.json] [--json OUT.json]

Example:
    python3 surface_props.py 1N8Z.cif C --atlas ../viewer-atlas/atlas.json \
        --json ../viewer-atlas/props.json

What it computes, per residue
-----------------------------
    rsasa        relative solvent-accessible surface, Shrake-Rupley, 1.4 A probe
    kd           Kyte-Doolittle hydrophobicity of the residue type, -4.5 to 4.5
    hydrophobic  exposed hydrophobic area: rsasa x kd, only where kd > 0, in A^2.
                 A buried leucine is not a hydrophobic patch; this says so.
    burial       antigen heavy atoms within 10 A of the residue centroid.
                 High means recessed, in a groove or pocket. Low means it sits
                 proud of the surface.

On "shape complementarity"
--------------------------
Sc, as Lawrence and Colman defined it, is a property of an *interface*: it
needs two surfaces facing each other. A lone antigen does not have one. What a
lone antigen has is local shape, and `burial` is the honest version of that:
a measure of how enclosed each patch is. Do not report it as Sc.

If an atlas of known binders is supplied, this also asks whether either
property predicts where binders actually land, by rank correlation against the
per-residue binder density. That is the part worth reading: a property that
does not correlate is not a design criterion, however appealing it looks.

Stdlib only.
"""
import sys, json, math, argparse
from collections import defaultdict

AA3 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
       'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
       'THR':'T','TRP':'W','TYR':'Y','VAL':'V','MSE':'M'}
# Kyte & Doolittle, J Mol Biol 157:105 (1982)
KD = {'A':1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C':2.5,'Q':-3.5,'E':-3.5,'G':-0.4,
      'H':-3.2,'I':4.5,'L':3.8,'K':-3.9,'M':1.9,'F':2.8,'P':-1.6,'S':-0.8,
      'T':-0.7,'W':-0.9,'Y':-1.3,'V':4.2}
# Tien et al. 2013, theoretical maximum ASA
MAXASA = {'A':129,'R':274,'N':195,'D':193,'C':167,'Q':225,'E':223,'G':104,'H':224,
          'I':197,'L':201,'K':236,'M':224,'F':240,'P':159,'S':155,'T':172,'W':285,
          'Y':263,'V':174}
RAD = {'C':1.70,'N':1.55,'O':1.52,'S':1.80,'SE':1.90}
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
        ic = r[ix['pdbx_PDB_ins_code']] if 'pdbx_PDB_ins_code' in ix else '?'
        out.append(dict(ch=r[ix['auth_asym_id']], rs=int(r[ix['auth_seq_id']]),
                        ic='' if ic in ('?', '.') else ic,
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


def spearman(xs, ys):
    """Rank correlation, with average ranks for ties."""
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    rx, ry = rank(xs), rank(ys)
    n = len(xs)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    dx = math.sqrt(sum((a - mx) ** 2 for a in rx))
    dy = math.sqrt(sum((b - my) ** 2 for b in ry))
    return num / (dx * dy) if dx and dy else 0.0


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('structure'); ap.add_argument('chain')
    ap.add_argument('--atlas')
    ap.add_argument('--json')
    args = ap.parse_args()

    atoms = [a for a in parse(args.structure) if a['ch'] == args.chain]
    if not atoms:
        sys.exit(f'chain {args.chain} not found')
    for a, v in zip(atoms, sasa(atoms)):
        a['sasa'] = v
    # Residues with an insertion code (52A beside 52) are not reported, as in atlas.py: keyed by
    # number alone they would be merged into the residue before them. Their atoms still count
    # towards the surface of their neighbours above.
    res = defaultdict(list)
    for a in atoms:
        if not a['ic']:
            res[a['rs']].append(a)
    ids = sorted(res)
    aa1 = {r: AA3.get(res[r][0]['aa'], 'X') for r in ids}

    def centroid(r):
        sc = [a for a in res[r] if a['at'] not in ('N', 'C', 'O')] or res[r]
        return tuple(sum(a['xyz'][k] for a in sc) / len(sc) for k in range(3))
    cen = {r: centroid(r) for r in ids}

    allxyz = [a['xyz'] for a in atoms]
    cell = 10.0
    grid = defaultdict(list)
    for n, (x, y, z) in enumerate(allxyz):
        grid[(int(x//cell), int(y//cell), int(z//cell))].append(n)

    props = {}
    for r in ids:
        a1 = aa1[r]
        area = sum(a['sasa'] for a in res[r])
        rs = area / MAXASA.get(a1, 200)
        kd = KD.get(a1, 0.0)
        p = cen[r]
        cx, cy, cz = (int(v//cell) for v in p); c = 0
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for n in grid.get((cx+dx, cy+dy, cz+dz), ()):
                        if math.dist(p, allxyz[n]) <= 10.0:
                            c += 1
        props[r] = {'aa': a1, 'area': round(area, 1), 'rsasa': round(rs, 3),
                    'kd': kd, 'hydrophobic': round(max(0.0, kd) * rs, 3),
                    'burial': c}

    print(f'{args.structure} chain {args.chain}: {len(ids)} residues')
    exposed = [r for r in ids if props[r]['rsasa'] >= 0.10]
    print(f'exposed (rSASA >= 0.10): {len(exposed)}')
    bl = sorted(props[r]['burial'] for r in ids)
    print(f'burial range {bl[0]}-{bl[-1]}, median {bl[len(bl)//2]}')
    hy = sorted((props[r]['hydrophobic'], r) for r in exposed)[-8:]
    print('most exposed hydrophobic residues: '
          + ', '.join(f'{props[r]["aa"]}{r}' for _, r in reversed(hy)))

    out = {'structure': args.structure, 'chain': args.chain,
           'scales': {'hydrophobicity': 'Kyte-Doolittle 1982',
                      'max_asa': 'Tien et al. 2013',
                      'burial': 'antigen heavy atoms within 10 A of the residue centroid'},
           'note': 'burial is local shape, not shape complementarity. Sc needs two '
                   'surfaces; a lone antigen does not have an interface.',
           'residues': {str(r): props[r] for r in ids}}

    if args.atlas:
        A = json.load(open(args.atlas))
        # atlas.py renamed this key to density_distinct. Reading only the
        # old name made every correlation read 0.000 and every residue
        # read "never bound", silently.
        raw = A.get('density_distinct') or A.get('density') or {}
        dens = {int(k): v for k, v in raw.items()}
        if not dens:
            print('  ! the atlas carries no density map; '
                  'the comparison below is empty')
        cd = {int(k): v for k, v in A.get('counts_distinct', {}).items()}
        xs_h = [props[r]['hydrophobic'] for r in exposed]
        xs_b = [props[r]['burial'] for r in exposed]
        xs_a = [props[r]['rsasa'] for r in exposed]
        ys_d = [dens.get(r, 0) for r in exposed]
        ys_c = [cd.get(r, 0) for r in exposed]
        cor = {
            'n_exposed': len(exposed),
            'hydrophobic_vs_density': round(spearman(xs_h, ys_d), 3),
            'burial_vs_density': round(spearman(xs_b, ys_d), 3),
            'rsasa_vs_density': round(spearman(xs_a, ys_d), 3),
            'hydrophobic_vs_binders': round(spearman(xs_h, ys_c), 3),
            'burial_vs_binders': round(spearman(xs_b, ys_c), 3),
            'rsasa_vs_binders': round(spearman(xs_a, ys_c), 3),
        }
        out['correlations'] = cor
        print('\nDoes surface property predict where binders land?')
        print('  Spearman rank correlation over the exposed residues, '
              f'n = {len(exposed)}')
        for k in ('hydrophobic', 'burial', 'rsasa'):
            print(f'    {k:<12} vs antibody-residue density: '
                  f'{cor[k + "_vs_density"]:+.3f}    vs distinct binders: '
                  f'{cor[k + "_vs_binders"]:+.3f}')

        # A correlation over mostly-zero data is hard to read. Compare the two
        # populations directly instead: residues some binder touches, against
        # residues none does.
        bound = [r for r in exposed if dens.get(r, 0) > 0]
        free = [r for r in exposed if dens.get(r, 0) == 0]
        def mean(v): return sum(v) / len(v) if v else 0.0
        comp = {'n_bound': len(bound), 'n_unbound': len(free)}
        print(f'\n  Bound ({len(bound)}) vs never bound ({len(free)}), exposed residues only:')
        print(f"    {'property':<14} {'bound':>9} {'unbound':>9}  reading")
        for key, label, how in (
                ('hydrophobic', 'hydrophobic', 'exposed hydrophobic area, rSASA x KD'),
                ('kd', 'Kyte-Doolittle', 'raw hydrophobicity of the residue type'),
                ('burial', 'burial', 'antigen atoms within 10 A'),
                ('rsasa', 'rSASA', 'how exposed it is')):
            b, f = mean([props[r][key] for r in bound]), mean([props[r][key] for r in free])
            comp[key] = {'bound': round(b, 3), 'unbound': round(f, 3)}
            print(f'    {label:<14} {b:>9.2f} {f:>9.2f}  {how}')
        out['bound_vs_unbound'] = comp

    if args.json:
        with open(args.json, 'w') as fh:
            json.dump(out, fh, separators=(',', ':'))
            fh.write('\n')
        print(f'\nwrote {args.json}')


if __name__ == '__main__':
    main()
