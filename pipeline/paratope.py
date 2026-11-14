#!/usr/bin/env python3
"""Measure an antibody's paratope and group it into loops.

Usage:
    python3 paratope.py STRUCTURE ANTIGEN_CHAIN AB_CHAINS [--cutoff 4.5]
                        [--naming B:H,A:L] [--gap 8] [--json OUT.json]

Example (trastuzumab Fab on HER2):
    python3 paratope.py 1N8Z.cif C AB --naming B:H,A:L --json paratope.json

What it does
------------
Finds every antibody residue with a heavy atom within the cutoff of the
antigen, counts the atom pairs, records which antigen residues it touches,
then splits each chain's contacts into loops wherever the sequence gap exceeds
`--gap`.

What it does not do
-------------------
It does not assign CDRs. It measures contacts and groups them. The loops it
finds usually line up with the CDRs, but Kabat, Chothia and IMGT draw CDR
boundaries differently, and none of them is what this script computed. If you
need assigned CDRs, run ANARCI and compare.

Stdlib only.
"""
import sys, json, math, argparse
from collections import defaultdict

AA3 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
       'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
       'THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
AROMATIC = set('YWF')


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
                        aa=r[ix['label_comp_id']],
                        xyz=(float(r[ix['Cartn_x']]), float(r[ix['Cartn_y']]),
                             float(r[ix['Cartn_z']]))))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('structure')
    ap.add_argument('antigen', help='antigen chain id')
    ap.add_argument('ab', help='antibody chain ids, e.g. AB')
    ap.add_argument('--cutoff', type=float, default=4.5)
    ap.add_argument('--gap', type=int, default=8,
                    help='a sequence gap larger than this starts a new loop')
    ap.add_argument('--naming', default='',
                    help='chain:prefix pairs, e.g. B:H,A:L')
    ap.add_argument('--json', help='write the result here')
    args = ap.parse_args()

    naming = {}
    for pair in filter(None, args.naming.split(',')):
        c, _, p = pair.partition(':')
        naming[c.strip()] = p.strip()

    atoms = parse(args.structure)
    ag = [a for a in atoms if a['ch'] == args.antigen]
    ab = [a for a in atoms if a['ch'] in list(args.ab)]
    if not ag or not ab:
        sys.exit(f'chains present: {sorted({a["ch"] for a in atoms})}')

    cell = args.cutoff
    grid = defaultdict(list)
    for n, a in enumerate(ag):
        x, y, z = a['xyz']
        grid[(int(x//cell), int(y//cell), int(z//cell))].append(n)

    res = {}
    for a in ab:
        x, y, z = a['xyz']
        cx, cy, cz = int(x//cell), int(y//cell), int(z//cell)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for n in grid.get((cx+dx, cy+dy, cz+dz), ()):
                        b = ag[n]
                        d = math.dist(a['xyz'], b['xyz'])
                        if d > args.cutoff:
                            continue
                        # The insertion code is part of the residue's identity. Kabat
                        # numbering has 52 and 52A. Keyed by number alone they merge.
                        key = (a['ch'], a['rs'], a['ic'])
                        e = res.setdefault(key, {'chain':a['ch'], 'resi':a['rs'],
                                                 'icode':a['ic'], 'label':f"{a['rs']}{a['ic']}",
                                                 'aa':AA3.get(a['aa'], 'X'),
                                                 'pairs':0, 'closest':99.0,
                                                 'touches':defaultdict(int)})
                        e['pairs'] += 1
                        e['closest'] = min(e['closest'], d)
                        e['touches'][b['rs']] += 1

    loops, by_chain = [], defaultdict(list)
    for key in sorted(res):
        by_chain[key[0]].append(res[key])
    for ch in sorted(by_chain):
        items = sorted(by_chain[ch], key=lambda e: (e['resi'], e['icode']))
        groups, cur = [], [items[0]]
        for e in items[1:]:
            if e['resi'] - cur[-1]['resi'] > args.gap:
                groups.append(cur); cur = []
            cur.append(e)
        groups.append(cur)
        prefix = naming.get(ch, ch)
        for n, g in enumerate(groups, 1):
            touches = sorted({t for e in g for t in e['touches']})
            loops.append({
                'id': f'{prefix}{n}', 'chain': ch,
                'res': sorted({e['resi'] for e in g}),          # numbers, for a viewer's selection
                'labels': [e['label'] for e in g],             # with insertion codes, for people
                'seq': ''.join(e['aa'] for e in g),
                'pairs': sum(e['pairs'] for e in g),
                'aromatics': [f"{e['aa']}{e['label']}" for e in g if e['aa'] in AROMATIC],
                'touches': touches,
                '_keys': [(e['chain'], e['resi'], e['icode']) for e in g]})

    residues = []
    for key in sorted(res, key=lambda k: -res[k]['pairs']):
        e = res[key]
        loop = next((l['id'] for l in loops if key in l['_keys']), None)
        residues.append({'chain':e['chain'], 'resi':e['resi'], 'icode':e['icode'],
                         'label':e['label'], 'aa':e['aa'],
                         'loop':loop, 'pairs':e['pairs'],
                         'closest':round(e['closest'], 2),
                         'touches':sorted(e['touches'])})

    for l in loops:
        l.pop('_keys')
    total = sum(r['pairs'] for r in residues)
    arom = sum(r['pairs'] for r in residues if r['aa'] in AROMATIC)
    out = {
        'structure': args.structure, 'antigen_chain': args.antigen,
        'antibody_chains': list(args.ab), 'cutoff': args.cutoff,
        'loop_gap': args.gap,
        'note': 'Loops are measured contact groups, not assigned CDRs. '
                'Kabat, Chothia and IMGT draw CDR boundaries differently.',
        'total_residues': len(residues), 'total_atom_pairs': total,
        'aromatic_residues': sum(1 for r in residues if r['aa'] in AROMATIC),
        'aromatic_atom_pairs': arom,
        'loops': loops, 'residues': residues,
    }

    print(f"{args.structure}: antigen {args.antigen} vs antibody "
          f"{'+'.join(args.ab)} at {args.cutoff} A")
    print(f"paratope: {len(residues)} residues, {total} atom pairs")
    print(f"aromatic (Y/W/F): {out['aromatic_residues']} residues carrying "
          f"{arom} pairs ({100*arom//total}%)\n")
    h = f"{'loop':>5} {'chain':>5} {'residues':<28} {'seq':<10} {'pairs':>6}  touches on the antigen"
    print(h); print('-' * min(len(h) + 26, 100))
    for l in sorted(loops, key=lambda l: -l['pairs']):
        print(f"{l['id']:>5} {l['chain']:>5} "
              f"{','.join(l['labels']):<28} {l['seq']:<10} {l['pairs']:>6}  "
              f"{','.join(str(t) for t in l['touches'])}")
    print(f"\n{'residue':>8} {'loop':>5} {'pairs':>6} {'closest':>8}  touches")
    for r in residues:
        print(f"{r['chain']}{r['aa']}{r['label']:<6} {r['loop']:>5} {r['pairs']:>6} "
              f"{r['closest']:>8}  {','.join(str(t) for t in r['touches'])}")

    if args.json:
        with open(args.json, 'w') as fh:
            json.dump(out, fh, indent=2)
            fh.write('\n')
        print(f"\nwrote {args.json}")


if __name__ == '__main__':
    main()
