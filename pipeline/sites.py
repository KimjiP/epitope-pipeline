#!/usr/bin/env python3
"""Find and rank candidate binder sites on an antigen in an mmCIF structure.

Usage:
    python3 sites.py STRUCTURE.cif ANTIGEN_CHAIN [REF_CHAINS|-] [LO] [HI]

Example (HER2 + trastuzumab Fab, domain IV only):
    python3 sites.py 1N8Z.cif C AB 490 607

What it does, and does not do
----------------------------
It computes real solvent-accessible surface area (Shrake-Rupley, 1.4 A probe)
for the antigen chain ALONE. It then groups exposed residues into patches of
10 A radius and scores each patch. If you name reference chains (an antibody
in the same file), it also measures how close each patch is to that antibody,
and how many residues the patch shares with that antibody's footprint.

Every threshold below is a choice, not a measurement:
    contact      = any heavy-atom pair within 4.5 A
    exposed      = relative SASA >= 0.25
    anchor       = hydrophobic (FWYLIVM) and relative SASA >= 0.30
    patch        = all exposed residues within 10 A of a seed side chain
    COMPETES     = a patch residue is within 12 A of the reference antibody
Change them and the answer changes. Report which values you used.

Stdlib only. No install step.
"""
import sys, math
from collections import defaultdict

AA3 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
       'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
       'THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
# Tien et al. 2013, theoretical maximum ASA per residue
MAXASA = {'A':129,'R':274,'N':195,'D':193,'C':167,'Q':225,'E':223,'G':104,'H':224,
          'I':197,'L':201,'K':236,'M':224,'F':240,'P':159,'S':155,'T':172,'W':285,
          'Y':263,'V':174}
RAD = {'C':1.70,'N':1.55,'O':1.52,'S':1.80,'SE':1.90}
HYDROPHOBIC = set('FWYLIVM')
PROBE, NPTS = 1.4, 96
CONTACT, EXP_MIN, ANCHOR_MIN, PATCH_R, COMPETE_A = 4.5, 0.25, 0.30, 10.0, 12.0


def parse(path):
    lines = open(path).read().splitlines()
    cols, rows, i = [], [], 0
    while i < len(lines):
        if lines[i].startswith('loop_'):
            j = i + 1; tags = []
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
        sys.exit("no _atom_site loop found - is this an mmCIF file?")
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
    """Shrake-Rupley per-atom SASA in A^2."""
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


def runs(sorted_ids):
    out, st, pv = [], sorted_ids[0], sorted_ids[0]
    for r in sorted_ids[1:]:
        if r - pv > 1:
            out.append((st, pv)); st = r
        pv = r
    out.append((st, pv))
    return ','.join(f'{a}-{b}' if a != b else str(a) for a, b in out)


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    path, AG = sys.argv[1], sys.argv[2]
    PARTNER = list(sys.argv[3]) if len(sys.argv) > 3 and sys.argv[3] != '-' else []
    LO = int(sys.argv[4]) if len(sys.argv) > 4 else -10**9
    HI = int(sys.argv[5]) if len(sys.argv) > 5 else 10**9

    atoms = parse(path)
    ag = [a for a in atoms if a['ch'] == AG]
    ab = [a for a in atoms if a['ch'] in PARTNER]
    if not ag:
        sys.exit(f"chain {AG} has no ATOM records. Chains present: "
                 f"{sorted({a['ch'] for a in atoms})}")

    res = defaultdict(list)
    for a in ag:
        res[a['rs']].append(a)
    resids = sorted(res)
    aa1 = {r: AA3.get(res[r][0]['aa'], 'X') for r in resids}

    for a, v in zip(ag, sasa(ag)):
        a['sasa'] = v
    rsa = {r: sum(a['sasa'] for a in res[r]) / MAXASA.get(aa1[r], 200) for r in resids}

    def centroid(r):
        sc = [a for a in res[r] if a['at'] not in ('N', 'C', 'O')] or res[r]
        return tuple(sum(a['xyz'][k] for a in sc) / len(sc) for k in range(3))
    cen = {r: centroid(r) for r in resids}

    footprint, ref_dist, pairs = set(), {}, defaultdict(int)
    if ab:
        abx = [a['xyz'] for a in ab]; cell = 6.0
        g2 = defaultdict(list)
        for n, (x, y, z) in enumerate(abx):
            g2[(int(x//cell), int(y//cell), int(z//cell))].append(n)

        def mind(p, span):
            cx, cy, cz = (int(v//cell) for v in p); best = 1e9
            for dx in range(-span, span+1):
                for dy in range(-span, span+1):
                    for dz in range(-span, span+1):
                        for n in g2.get((cx+dx, cy+dy, cz+dz), ()):
                            d = math.dist(p, abx[n])
                            if d < best: best = d
            return best
        for r in resids:
            for a in res[r]:
                if mind(a['xyz'], 1) <= CONTACT:
                    footprint.add(r); pairs[r] += 1
            ref_dist[r] = mind(cen[r], 6)

    glyc = [r for r in resids if aa1[r] == 'N'
            and aa1.get(r+1) != 'P' and aa1.get(r+2) in ('S', 'T')]
    cys = [r for r in resids if aa1[r] == 'C']
    gaps = [(a+1, b-1) for a, b in zip(resids, resids[1:]) if b - a > 1]

    exposed = {r for r in resids if rsa[r] >= EXP_MIN}
    pool = [r for r in resids if LO <= r <= HI and r in exposed]
    patches = []
    for s in pool:
        mem = [r for r in pool if math.dist(cen[s], cen[r]) <= PATCH_R]
        if len(mem) < 5:
            continue
        anchors = [r for r in mem if aa1[r] in HYDROPHOBIC and rsa[r] >= ANCHOR_MIN]
        patches.append(dict(
            seed=s, mem=mem, n=len(mem), anchors=anchors,
            rsa=sum(rsa[r] for r in mem) / len(mem),
            area=sum(sum(a['sasa'] for a in res[r]) for r in mem),
            gl=min((math.dist(cen[r], cen[g]) for r in mem for g in glyc), default=999),
            cy=min((math.dist(cen[r], cen[c]) for r in mem for c in cys), default=999),
            ref=(min(ref_dist[r] for r in mem) if ab else None),
            ov=len(set(mem) & footprint)))

    def score(p):
        s = p['area'] / 100 + 3 * len(p['anchors']) + 6 * p['rsa']
        if p['gl'] < 12: s -= 8
        if p['cy'] < 10: s += 3
        return s
    patches.sort(key=score, reverse=True)
    chosen = []
    for p in patches:
        if all(len(set(p['mem']) & set(q['mem'])) < 0.4 * min(p['n'], q['n'])
               for q in chosen):
            chosen.append(p)
        if len(chosen) >= 10:
            break

    print(f"file            {path}")
    print(f"antigen chain   {AG}: {len(resids)} modelled residues, {resids[0]}-{resids[-1]}")
    print(f"unmodelled      {gaps or 'none'}   <-- disordered or absent; a design "
          f"cannot use these")
    print(f"N-glyc sequons  {[f'N{r}' for r in glyc] or 'none'}   (sequence motif N-X-S/T; "
          f"not proof of a glycan)")
    if ab:
        print(f"reference       chains {'+'.join(PARTNER)}")
        print(f"its footprint   {runs(sorted(footprint))}   [{len(footprint)} residues]")
    print(f"patch window    {max(LO, resids[0])}-{min(HI, resids[-1])}")
    print(f"thresholds      contact {CONTACT} A | exposed rSASA >= {EXP_MIN} | "
          f"anchor rSASA >= {ANCHOR_MIN} | patch r {PATCH_R} A | competes <= {COMPETE_A} A\n")

    h = (f"{'#':>2} {'seed':>6} {'n':>3} {'area A2':>8} {'rSASA':>6} {'anch':>5} "
         f"{'glyc_A':>7} {'ref_A':>6} {'ovlp':>5}  verdict")
    print(h); print('-' * len(h))
    for n, p in enumerate(chosen, 1):
        v = []
        if p['ref'] is not None:
            v.append('COMPETES' if p['ref'] <= COMPETE_A
                     else ('rim' if p['ref'] <= 20 else 'free site'))
        if p['gl'] < 12: v.append('glycan risk')
        if len(p['anchors']) >= 3: v.append('anchors')
        if p['cy'] < 10: v.append('rigid')
        rd = '>30' if (p['ref'] is None or p['ref'] > 30) else f"{p['ref']:.1f}"
        print(f"{n:>2} {aa1[p['seed']]}{p['seed']:<5} {p['n']:>3} {p['area']:>8.0f} "
              f"{p['rsa']:>6.2f} {len(p['anchors']):>5} {p['gl']:>7.1f} {rd:>6} "
              f"{p['ov']:>5}  {'; '.join(v)}")

    print()
    for n, p in enumerate(chosen, 1):
        print(f"[{n}] seed {aa1[p['seed']]}{p['seed']}  segments {runs(sorted(p['mem']))}")
        print(f"     anchors      {','.join(f'{aa1[r]}{r}' for r in p['anchors']) or '(none)'}")
        print(f"     in footprint {','.join(f'{aa1[r]}{r}' for r in sorted(set(p['mem']) & footprint)) or '-'}")

    keep = sorted(footprint | {r for p in chosen for r in p['anchors']})
    if keep:
        print(f"\nper-residue detail (footprint residues and all anchors)")
        print(f"{'res':>7} {'rSASA':>6} {'area A2':>8} {'pairs':>6} {'ref_A':>6}  flags")
        for r in keep:
            fl = []
            if r in footprint: fl.append('footprint')
            if aa1[r] in HYDROPHOBIC and rsa[r] >= ANCHOR_MIN: fl.append('anchor')
            if r in glyc: fl.append('N-glyc sequon')
            if aa1[r] == 'C': fl.append('cysteine')
            rd = f"{ref_dist[r]:.1f}" if ab and ref_dist.get(r, 1e9) <= 30 else '-'
            print(f"{aa1[r]}{r:<6} {rsa[r]:>6.2f} "
                  f"{sum(a['sasa'] for a in res[r]):>8.0f} {pairs.get(r, 0):>6} {rd:>6}  "
                  f"{', '.join(fl)}")


if __name__ == '__main__':
    main()
