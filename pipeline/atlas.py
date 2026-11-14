#!/usr/bin/env python3
"""Build an epitope atlas: every distinct binder of one antigen, in one frame.

Usage:
    python3 atlas.py REF.cif REF_CHAIN REFERENCE.fasta CACHE_DIR
        [--exclude FASTA] [--entities JSON] [--name NAME=PDB ...]
        [--cutoff 4.5] [--near 5.0] [--identity 0.95] [--min-cover 50]
        [--domains "I=1-195,II=196-319,..."] [--out DIR] [--selftest]

Example:
    python3 atlas.py 1N8Z.cif C P04626.fasta ../pdb-cache --out ../atlas-out \
        --exclude erbb_others.fasta --entities ../pdb-cache/entities.json \
        --name trastuzumab=1N8Z --name pertuzumab=1S78

What it does
------------
1. Reads every .cif in CACHE_DIR.
2. Finds the antigen chain in each, by aligning chains to the reference
   sequence. This also recovers each entry's own numbering offset, which no
   file records.
3. Finds the other chains that touch the antigen. Chains of relatives of the
   antigen (--exclude) are receptor partners, not binders, and are dropped.
4. Splits the binder chains into UNITS. A unit is a set of chains that touch
   each other: one Fab, one nanobody, one scFv. A complex with two Fabs gives
   two units. Counting entries instead of units counts one antibody once for
   every entry it appears in, and counts a two-antibody complex as one binder.
5. Superposes the antigen onto the reference (Horn's quaternion method) and
   moves every unit into the reference frame.
6. Groups units into DISTINCT BINDERS by sequence identity (Needleman-Wunsch).
   Antibody chains are compared over the variable domain only, because two
   different antibodies on one framework share most of their sequence.
7. Counts, in the reference numbering, for every antigen residue:
      counts_units     how many units touch it
      counts_distinct  how many distinct binders touch it
      density_distinct antibody residues within --near of it, summed over
                       distinct binders (the most any one structure shows)
8. With --name, builds a footprint for a named binder from ALL of its
   structures. One structure can miss a loop that another one models.

Everything is in the reference's numbering, so entries that use different
numbering conventions can be compared directly.

Stdlib only.
"""
import sys, os, json, math, glob, argparse, random, string, datetime
from collections import defaultdict

# one single-character chain id per unit, so a viewer can select each binder
CODES = string.ascii_uppercase + string.ascii_lowercase + string.digits

AA3 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
       'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
       'THR':'T','TRP':'W','TYR':'Y','VAL':'V','MSE':'M'}


# ------------------------------------------------------------------ parse --

def parse(path):
    """ATOM records of model 1 as (chain, resi, icode, comp, atom, xyz)."""
    try:
        lines = open(path, errors='replace').read().splitlines()
    except OSError:
        return []
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
        return []
    ix = {c: n for n, c in enumerate(cols)}
    need = ('group_PDB', 'auth_asym_id', 'auth_seq_id', 'label_comp_id',
            'label_atom_id', 'type_symbol', 'Cartn_x', 'Cartn_y', 'Cartn_z')
    if any(n not in ix for n in need):
        return []
    model1 = None
    out = []
    for r in rows:
        if r[ix['group_PDB']] != 'ATOM' or r[ix['type_symbol']] == 'H':
            continue
        if 'pdbx_PDB_model_num' in ix:
            m = r[ix['pdbx_PDB_model_num']]
            if model1 is None:
                model1 = m
            elif m != model1:
                continue
        alt = r[ix['label_alt_id']] if 'label_alt_id' in ix else '.'
        if alt not in ('.', '?', '', 'A'):
            continue
        try:
            rs = int(r[ix['auth_seq_id']])
        except ValueError:
            continue
        ic = r[ix['pdbx_PDB_ins_code']] if 'pdbx_PDB_ins_code' in ix else '?'
        ic = '' if ic in ('?', '.') else ic
        out.append((r[ix['auth_asym_id']], rs, ic, r[ix['label_comp_id']],
                    r[ix['label_atom_id']],
                    (float(r[ix['Cartn_x']]), float(r[ix['Cartn_y']]),
                     float(r[ix['Cartn_z']]))))
    return out


def chain_residues(atoms):
    """{chain: {(resi, icode): (comp, {atom: xyz})}}

    The key keeps the insertion code. Kabat-numbered antibodies carry residues
    such as 52 and 52A. Keyed by number alone, the two collapse into one
    residue and one of them loses its atoms."""
    ch = defaultdict(dict)
    for c, rs, ic, comp, at, xyz in atoms:
        d = ch[c].setdefault((rs, ic), (comp, {}))
        d[1][at] = xyz
    return ch


def intkeyed(res):
    """Residues without an insertion code, keyed by number. For the antigen."""
    return {k[0]: v for k, v in res.items() if k[1] == ''}


def seq_of(res):
    return ''.join(AA3.get(res[k][0], 'X') for k in sorted(res))


# -------------------------------------------------------------- alignment --

def find_offset(res, ref):
    """Offset o such that reference index (1-based) = resi + o, or None.
    `res` is {int resi: (comp, atoms)}."""
    ids = sorted(res)
    seq = {r: AA3.get(res[r][0], 'X') for r in ids}
    windows = []
    for n in range(len(ids) - 11):
        block = ids[n:n + 12]
        if block[-1] - block[0] == 11:
            windows.append(block)
        if len(windows) >= 12:
            break
    best = (0.0, None)
    tried = set()
    for block in windows:
        probe = ''.join(seq[r] for r in block)
        if 'X' in probe:
            continue
        start = 0
        while True:
            hit = ref.find(probe, start)
            if hit < 0:
                break
            o = (hit + 1) - block[0]
            start = hit + 1
            if o in tried:
                continue
            tried.add(o)
            ok = tot = 0
            for r in ids:
                p = r + o
                if 1 <= p <= len(ref):
                    tot += 1
                    if seq[r] == ref[p - 1]:
                        ok += 1
            if tot >= 20:
                score = ok / tot
                if score > best[0]:
                    best = (score, o)
    return best[1] if best[0] >= 0.80 else None


def nw_identity(a, b, match=2, mismatch=-1, gap=-2):
    """Global alignment. Identity = matched positions / length of the longer
    sequence, so a truncated copy scores lower than a full one."""
    n, m = len(a), len(b)
    if n == 0 or m == 0:
        return 0.0
    ptr = [[0] * (m + 1) for _ in range(n + 1)]
    prev = [j * gap for j in range(m + 1)]
    for j in range(1, m + 1):
        ptr[0][j] = 2
    for i in range(1, n + 1):
        cur = [i * gap] + [0] * m
        ptr[i][0] = 1
        ai = a[i - 1]
        row = ptr[i]
        for j in range(1, m + 1):
            s = prev[j - 1] + (match if ai == b[j - 1] else mismatch)
            u = prev[j] + gap
            l = cur[j - 1] + gap
            if s >= u and s >= l:
                cur[j] = s; row[j] = 0
            elif u >= l:
                cur[j] = u; row[j] = 1
            else:
                cur[j] = l; row[j] = 2
        prev = cur
    i, j, hits = n, m, 0
    while i > 0 or j > 0:
        p = ptr[i][j]
        if i > 0 and j > 0 and p == 0:
            if a[i - 1] == b[j - 1]:
                hits += 1
            i -= 1; j -= 1
        elif i > 0 and (p == 1 or j == 0):
            i -= 1
        else:
            j -= 1
    return hits / max(n, m)


# ---------------------------------------------------------- superposition --

def jacobi(a, n):
    a = [row[:] for row in a]
    v = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    for _ in range(100):
        off = sum(a[i][j] ** 2 for i in range(n) for j in range(n) if i != j)
        if off < 1e-20:
            break
        p, q = max(((i, j) for i in range(n) for j in range(n) if i < j),
                   key=lambda t: abs(a[t[0]][t[1]]))
        if abs(a[p][q]) < 1e-20:
            break
        theta = (a[q][q] - a[p][p]) / (2 * a[p][q])
        t = (1 if theta >= 0 else -1) / (abs(theta) + math.sqrt(theta * theta + 1))
        c = 1 / math.sqrt(t * t + 1); s = t * c
        for k in range(n):
            akp, akq = a[k][p], a[k][q]
            a[k][p] = c * akp - s * akq
            a[k][q] = s * akp + c * akq
        for k in range(n):
            apk, aqk = a[p][k], a[q][k]
            a[p][k] = c * apk - s * aqk
            a[q][k] = s * apk + c * aqk
        for k in range(n):
            vkp, vkq = v[k][p], v[k][q]
            v[k][p] = c * vkp - s * vkq
            v[k][q] = s * vkp + c * vkq
    vals = [a[i][i] for i in range(n)]
    vecs = [[v[k][i] for k in range(n)] for i in range(n)]
    return vals, vecs


def superpose(mov, ref):
    """Rotation+translation putting `mov` onto `ref`. Horn's quaternion method.
    Returns (R, t, rmsd); apply as R*(x) + t."""
    n = len(mov)
    mc = [sum(p[i] for p in mov) / n for i in range(3)]
    rc = [sum(p[i] for p in ref) / n for i in range(3)]
    S = [[0.0] * 3 for _ in range(3)]
    for a, b in zip(mov, ref):
        da = (a[0] - mc[0], a[1] - mc[1], a[2] - mc[2])
        db = (b[0] - rc[0], b[1] - rc[1], b[2] - rc[2])
        for i in range(3):
            for j in range(3):
                S[i][j] += da[i] * db[j]
    xx, xy, xz = S[0]; yx, yy, yz = S[1]; zx, zy, zz = S[2]
    N = [[xx + yy + zz, yz - zy,        zx - xz,        xy - yx],
         [yz - zy,      xx - yy - zz,   xy + yx,        zx + xz],
         [zx - xz,      xy + yx,       -xx + yy - zz,   yz + zy],
         [xy - yx,      zx + xz,        yz + zy,       -xx - yy + zz]]
    vals, vecs = jacobi(N, 4)
    q = vecs[max(range(4), key=lambda i: vals[i])]
    w, x, y, z = q
    nrm = math.sqrt(w*w + x*x + y*y + z*z) or 1.0
    w, x, y, z = w/nrm, x/nrm, y/nrm, z/nrm
    R = [[w*w + x*x - y*y - z*z, 2*(x*y - w*z),         2*(x*z + w*y)],
         [2*(x*y + w*z),         w*w - x*x + y*y - z*z, 2*(y*z - w*x)],
         [2*(x*z - w*y),         2*(y*z + w*x),         w*w - x*x - y*y + z*z]]
    t = [rc[i] - sum(R[i][j] * mc[j] for j in range(3)) for i in range(3)]
    se = 0.0
    for a, b in zip(mov, ref):
        p = [sum(R[i][j] * a[j] for j in range(3)) + t[i] for i in range(3)]
        se += sum((p[i] - b[i]) ** 2 for i in range(3))
    return R, t, math.sqrt(se / n)


def apply_rt(R, t, p):
    return (sum(R[0][j] * p[j] for j in range(3)) + t[0],
            sum(R[1][j] * p[j] for j in range(3)) + t[1],
            sum(R[2][j] * p[j] for j in range(3)) + t[2])


def selftest():
    ok = True
    random.seed(7)
    pts = [(random.uniform(-20, 20), random.uniform(-20, 20), random.uniform(-20, 20))
           for _ in range(40)]
    R0, t0, rm = superpose(pts, pts)
    ident = max(abs(R0[i][j] - (1.0 if i == j else 0.0)) for i in range(3) for j in range(3))
    print(f"  superposition, identity case : rmsd {rm:.1e}, max deviation from I {ident:.1e}")
    ang = 0.7
    Rz = [[math.cos(ang), -math.sin(ang), 0], [math.sin(ang), math.cos(ang), 0], [0, 0, 1]]
    shifted = [apply_rt(Rz, [5.0, -3.0, 2.0], p) for p in pts]
    R1, t1, rm1 = superpose(shifted, pts)
    print(f"  superposition, known rotation: rmsd {rm1:.1e}")
    ok &= rm < 1e-6 and ident < 1e-6 and rm1 < 1e-6

    a = 'MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVV'
    L = len(a)
    b = a[:10] + ('W' if a[10] != 'W' else 'Y') + a[11:]
    c = a[:20] + a[25:]
    d = ''.join(random.choice('ACDEFGHIKLMNPQRSTVWY') for _ in range(L))
    i1, i2, i3, i4 = nw_identity(a, a), nw_identity(a, b), nw_identity(a, c), nw_identity(a, d)
    print(f"  identity, same sequence      : {i1:.3f}   (expect 1.000)")
    print(f"  identity, one substitution   : {i2:.3f}   (expect {(L-1)/L:.3f})")
    print(f"  identity, 5-residue deletion : {i3:.3f}   (expect {(L-5)/L:.3f})")
    print(f"  identity, unrelated sequence : {i4:.3f}   (expect well below 0.5)")
    ok &= abs(i1 - 1) < 1e-12 and abs(i2 - (L-1)/L) < 1e-12
    ok &= abs(i3 - (L-5)/L) < 1e-12 and i4 < 0.5
    return ok


# ------------------------------------------------------------- geometry ----

def build_grid(points, cell):
    g = defaultdict(list)
    for n, p in enumerate(points):
        g[(int(p[0] // cell), int(p[1] // cell), int(p[2] // cell))].append(n)
    return g


def neighbours(g, p, cell):
    cx, cy, cz = int(p[0] // cell), int(p[1] // cell), int(p[2] // cell)
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for dz in (-1, 0, 1):
                yield from g.get((cx + dx, cy + dy, cz + dz), ())


def interface(ag_atoms, unit_atoms, cut, near_cut):
    """ag_atoms: [(antigen residue, xyz)]. unit_atoms: [(unit residue id, xyz)].
    Returns pairs{antigen residue: heavy-atom pairs within cut} and
    near{antigen residue: set of unit residues within near_cut}."""
    cell = max(cut, near_cut)
    g = build_grid([p for _, p in unit_atoms], cell)
    pairs, near = defaultdict(int), defaultdict(set)
    for r, p in ag_atoms:
        for n in neighbours(g, p, cell):
            d = math.dist(p, unit_atoms[n][1])
            if d <= near_cut:
                near[r].add(unit_atoms[n][0])
                if d <= cut:
                    pairs[r] += 1
    return pairs, near


def chain_pairs(pa, pb, cut, stop_at):
    """Heavy-atom pairs within cut between two point lists, stopping early."""
    if len(pa) > len(pb):
        pa, pb = pb, pa
    g = build_grid(pb, cut)
    n = 0
    for p in pa:
        for k in neighbours(g, p, cut):
            if math.dist(p, pb[k]) <= cut:
                n += 1
                if n >= stop_at:
                    return n
    return n


# ----------------------------------------------------------- classification --

def short_name(descs):
    out = []
    for d in descs:
        s = d
        for junk in ('Variable Heavy chain of ', 'Variable Light chain of ',
                     ' Fragment-Light Chain', ' Fragment-Heavy Chain',
                     ' FAB LIGHT CHAIN', ' FAB HEAVY CHAIN', ' FAB light chain', ' FAB heavy chain',
                     ' - heavy chain', ' - light chain', ' heavy chain', ' light chain',
                     ' Heavy Chain', ' Light Chain', ' HEAVY CHAIN', ' LIGHT CHAIN',
                     ' H chain', ' L chain', ' H Chain', ' L Chain', ' HC', ' LC', ' VH-CH1'):
            s = s.replace(junk, '')
        s = s.strip(' -')
        if s and s not in out:
            out.append(s)
    return ' + '.join(out) or 'unnamed'


def classify(descs, title, obs_lens):
    n = len(obs_lens)
    d = ' '.join(descs).lower()
    t = (title or '').lower()
    if 'ankyrin' in d or 'darpin' in d:
        return 'DARPin'
    if 'repebody' in d:
        return 'Repebody'
    if 'nanobody' in d or 'vhh' in d:
        return 'VHH'
    if 'immunoglobulin g-binding protein a' in d or 'affibody' in d:
        return 'Affibody'
    if 'ig gamma' in d and n == 1:
        return 'Fcab (engineered Fc)'
    if 'scfv' in d or 'single-chain' in d:
        return 'scFv'
    if n == 2 and any(w in d for w in ('fab', 'heavy', 'light', 'antibody', 'igg', 'vh', 'vl')):
        return 'Fab'
    if n == 2 and all(170 <= L <= 270 for L in obs_lens):
        return 'Fab'
    if n == 1 and ('heavy chain' in d or 'light chain' in d):
        return 'antibody fragment'
    if n == 1:
        L = obs_lens[0]
        if 'scfv' in t or 'single-chain' in t:
            return 'scFv'
        if L <= 80:
            return 'small scaffold'
        if L <= 150:
            return 'single domain (VHH-sized)'
        if 200 <= L <= 300:
            return 'scFv'
    return 'other'


ANTIBODY_LIKE = ('Fab', 'scFv', 'antibody fragment')
VARIABLE_DOMAIN = 125     # residues compared for antibody chains


def comparable_seq(seq, fmt):
    """The part of a chain that carries the binding loops.
    Antibody chains: the variable domain. Fcab: the CH3 domain, where the loops
    were engineered (the CH2 domain is identical in every Fcab, and over the whole Fc two
    different Fcabs look 97% the same). Other scaffolds: the whole chain."""
    if fmt in ANTIBODY_LIKE:
        return seq[:VARIABLE_DOMAIN]
    if fmt.startswith('Fcab'):
        return seq[-110:]
    return seq


def unit_identity(ua, ub):
    """Mean identity of the best chain pairing, or 0 when the units differ in
    format or chain count."""
    if ua['format'] != ub['format'] or len(ua['seqs']) != len(ub['seqs']):
        return 0.0
    A = [comparable_seq(s, ua['format']) for s in ua['seqs']]
    B = [comparable_seq(s, ub['format']) for s in ub['seqs']]
    for a, b in zip(sorted(A, key=len), sorted(B, key=len)):
        if min(len(a), len(b)) < 0.6 * max(len(a), len(b)):
            return 0.0
    if len(A) == 1:
        return nw_identity(A[0], B[0])
    if len(A) == 2:
        straight = (nw_identity(A[0], B[0]) + nw_identity(A[1], B[1])) / 2
        crossed = (nw_identity(A[0], B[1]) + nw_identity(A[1], B[0])) / 2
        return max(straight, crossed)
    return sum(nw_identity(a, b) for a, b in zip(A, B)) / len(A)


def load_entities(path):
    info = {}
    for e in json.load(open(path))['data']['entries']:
        ri = e.get('rcsb_entry_info') or {}
        res = ri.get('resolution_combined') or [None]
        desc = {}
        for pe in e.get('polymer_entities') or []:
            name = (pe.get('rcsb_polymer_entity') or {}).get('pdbx_description') or ''
            for c in pe['rcsb_polymer_entity_container_identifiers']['auth_asym_ids']:
                desc[c] = name
        info[e['rcsb_id']] = {'title': (e.get('struct') or {}).get('title', ''),
                              'resolution': res[0], 'method': ri.get('experimental_method'),
                              'desc': desc}
    return info


# ------------------------------------------------------------------- main --

def runs_text(ids):
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


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('ref', nargs='?'); ap.add_argument('ref_chain', nargs='?')
    ap.add_argument('fasta', nargs='?'); ap.add_argument('cache', nargs='?')
    ap.add_argument('--cutoff', type=float, default=4.5, help='contact distance, A')
    ap.add_argument('--near', type=float, default=5.0, help='density distance, A')
    ap.add_argument('--identity', type=float, default=0.95,
                    help='sequence identity at which two units are one binder')
    ap.add_argument('--min-cover', type=int, default=50)
    ap.add_argument('--partner-pairs', type=int, default=10,
                    help='heavy-atom pairs with the antigen for a chain to count as a binder')
    ap.add_argument('--link-pairs', type=int, default=60,
                    help='heavy-atom pairs between two chains for them to be one unit '
                         '(a Fab heavy and light chain share hundreds)')
    ap.add_argument('--copy-ratio', type=float, default=0.5,
                    help='within one entry, an identical copy with fewer than this fraction '
                         'of the best copy\'s contacts is a crystal-packing contact, dropped')
    ap.add_argument('--footprint-identity', type=float, default=0.98,
                    help='a named binder\'s footprint uses only members at least this '
                         'identical to it; engineered variants are listed apart')
    ap.add_argument('--out', default='atlas-out')
    ap.add_argument('--only', metavar='IDS',
                    help='comma-separated PDB ids; use only these entries from the\n                          cache directory, whatever else is in it')
    ap.add_argument('--selftest', action='store_true')
    ap.add_argument('--exclude', metavar='FASTA',
                    help='sequences whose chains are NOT binders: relatives of the '
                         'antigen that form receptor complexes with it')
    ap.add_argument('--entities', metavar='JSON',
                    help='RCSB GraphQL result with the real chain names and resolutions')
    ap.add_argument('--sites', metavar='JSON',
                    help='candidate sites {name: [residues]}, copied into the output')
    ap.add_argument('--name', action='append', default=[], metavar='NAME=PDB',
                    help='name the binder found in this entry; builds its footprint')
    ap.add_argument('--domains', default='I=1-195,II=196-319,III=320-488,IV=489-607',
                    help='antigen domains for the reference numbering, or "" for none')
    args = ap.parse_args()

    if args.selftest:
        print("self-test")
        ok = selftest()
        print("  PASS" if ok else "  FAIL")
        if not ok:
            return 2
        print()
        if args.ref is None:        # --selftest on its own is the whole job
            return 0

    domains = []
    for part in filter(None, args.domains.split(',')):
        nm, _, rg = part.partition('=')
        lo, _, hi = rg.partition('-')
        domains.append((nm.strip(), int(lo), int(hi)))

    def domain_of(r):
        for nm, lo, hi in domains:
            if lo <= r <= hi:
                return nm
        return None

    excl = []
    if args.exclude:
        name, buf = None, []
        for line in open(args.exclude):
            if line.startswith('>'):
                if name:
                    excl.append((name, ''.join(buf)))
                name, buf = (line[1:].split('|')[2].split()[0] if '|' in line
                             else line[1:40].strip()), []
            else:
                buf.append(line.strip())
        if name:
            excl.append((name, ''.join(buf)))
        print(f"excluding chains that match: {', '.join(n for n, _ in excl)}")
    ents = load_entities(args.entities) if args.entities else {}
    if ents:
        print(f"chain names and resolutions for {len(ents)} entries from {args.entities}")
    named_in = {}
    for spec in args.name:
        nm, _, pdb = spec.partition('=')
        named_in[nm.strip()] = pdb.strip().upper()

    if not all((args.ref, args.ref_chain, args.fasta, args.cache)):
        sys.exit('need REF.cif REF_CHAIN REFERENCE.fasta CACHE_DIR (or --selftest alone)')
    ref_seq = ''.join(l.strip() for l in open(args.fasta) if not l.startswith('>'))
    ref_ch = intkeyed(chain_residues(parse(args.ref))[args.ref_chain])
    ref_off = find_offset(ref_ch, ref_seq)
    if ref_off is None:
        sys.exit('could not align the reference chain')
    ref_ca = {rs: ats['CA'] for rs, (comp, ats) in ref_ch.items() if 'CA' in ats}
    print(f"reference {os.path.basename(args.ref)} chain {args.ref_chain}: "
          f"{len(ref_ca)} CA atoms, offset {ref_off:+d}\n")

    os.makedirs(args.out, exist_ok=True)
    entries, units = [], []
    trace_lines, serial = [], 0
    modelled_by_entry = {}

    # Which entries this run may see. Without --only the whole cache dir is
    # used, so a cache shared between antigens silently lets another antigen's
    # structures into the atlas. run_pipeline.py always passes the exact list.
    paths = sorted(glob.glob(os.path.join(args.cache, '*.cif')))
    if args.only:
        want = {x.strip().upper() for x in args.only.split(',') if x.strip()}
        paths = [p for p in paths
                 if os.path.basename(p)[:-4].upper() in want]
        missing = want - {os.path.basename(p)[:-4].upper() for p in paths}
        if missing:
            sys.exit('--only names entries that are not in the cache: '
                     + ', '.join(sorted(missing)))
    for path in paths:
        pdb = os.path.basename(path).split('.')[0].upper()
        atoms = parse(path)
        if not atoms:
            entries.append({'pdb': pdb, 'status': 'unreadable'}); continue
        ch = chain_residues(atoms)
        meta = ents.get(pdb, {})

        ag, offs = [], {}
        for c, res in ch.items():
            ir = intkeyed(res)
            if len(ir) < 30:
                continue
            o = find_offset(ir, ref_seq)
            if o is None:
                continue
            cover = sum(1 for r in ir if (r + o - ref_off) in ref_ca)
            if cover >= args.min_cover:
                ag.append(c); offs[c] = o
        if not ag:
            entries.append({'pdb': pdb, 'status': 'no chain covers the reference region'})
            continue

        c = max(ag, key=lambda x: sum(1 for r in intkeyed(ch[x]) if (r + offs[x] - ref_off) in ref_ca))
        o = offs[c]
        shift = o - ref_off          # entry numbering + shift = reference numbering
        agres = intkeyed(ch[c])
        ag_atoms = [(r + shift, p) for r, (comp, ats) in agres.items() for p in ats.values()]
        modelled = sorted(r + shift for r in agres)

        # candidate chains: not the antigen, not a relative of the antigen
        cand, dropped = {}, []
        for c2, res2 in ch.items():
            if c2 == c or c2 in ag or len(res2) < 25:
                continue
            rel = next((n for n, sq in excl if find_offset(intkeyed(res2), sq) is not None), None)
            if rel:
                dropped.append(f'{c2}:{rel}')
                continue
            cand[c2] = [((c2, k), p) for k, (comp, ats) in res2.items() for p in ats.values()]
        partners = {}
        for c2, atoms2 in cand.items():
            pr, _ = interface(ag_atoms, atoms2, args.cutoff, args.cutoff)
            if sum(pr.values()) >= args.partner_pairs:
                partners[c2] = atoms2
        if not partners:
            entries.append({'pdb': pdb,
                            'status': ('only receptor partners: ' + ','.join(dropped))
                                      if dropped else 'antigen present, no contacting partner',
                            'antigen_chain': c, 'offset': o, 'dropped_partners': dropped})
            continue

        # units: chains that touch each other. A light chain that never touches the
        # antigen still belongs to its Fab, so link through every candidate chain.
        parent = {n: n for n in cand}
        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]; x = parent[x]
            return x
        pts = {n: [p for _, p in cand[n]] for n in cand}
        names = sorted(cand)
        for i, a in enumerate(names):
            for b in names[i + 1:]:
                if a in partners or b in partners:
                    if chain_pairs(pts[a], pts[b], args.cutoff, args.link_pairs) >= args.link_pairs:
                        parent[find(a)] = find(b)
        groups = defaultdict(list)
        for n in names:
            groups[find(n)].append(n)
        groups = {k: v for k, v in groups.items() if any(m in partners for m in v)}

        sh_pairs = [(ats['CA'], ref_ca[r + shift]) for r, (comp, ats) in agres.items()
                    if (r + shift) in ref_ca and 'CA' in ats]
        if len(sh_pairs) < 30:
            entries.append({'pdb': pdb, 'status': 'too few shared CA atoms to superpose'})
            continue
        R, t, rmsd = superpose([p[0] for p in sh_pairs], [p[1] for p in sh_pairs])
        modelled_by_entry[pdb] = set(modelled)

        proto = []
        for members in sorted(groups.values(), key=lambda g: min(g)):
            members = sorted(members, key=lambda x: -len(ch[x]))
            uatoms = [a for m in members for a in cand[m]]
            pairs, near = interface(ag_atoms, uatoms, args.cutoff, args.near)
            descs = [meta.get('desc', {}).get(m, '') for m in members]
            obs = [len(ch[m]) for m in members]
            proto.append({'members': members, 'pairs': pairs, 'near': near, 'descs': descs,
                          'obs': obs, 'format': classify(descs, meta.get('title'), obs),
                          'seqs': [seq_of(ch[m]) for m in members],
                          'n_pairs': sum(pairs.values())})

        # identical copies in one entry: only the copy that really binds is a binder
        # Identical copies in one entry: only the copy that really binds is a
        # binder. Index by position and re-read the list each time. Testing the
        # loop variable instead let an entry with three or more identical
        # copies drop the same one twice, which raised ValueError on an EGFR
        # entry carrying four copies of nanobody 7D12.
        gone = []
        for i in range(len(proto)):
            for j in range(i + 1, len(proto)):
                a, b = proto[i], proto[j]
                if a is None or b is None or unit_identity(a, b) < 0.95:
                    continue
                wi, si = (i, j) if a['n_pairs'] < b['n_pairs'] else (j, i)
                weak, strong = proto[wi], proto[si]
                if weak['n_pairs'] < args.copy_ratio * strong['n_pairs']:
                    gone.append({'chains': weak['members'], 'pairs': weak['n_pairs'],
                                 'strong_pairs': strong['n_pairs']})
                    proto[wi] = None
                    if wi == i:
                        break
        proto = [u for u in proto if u is not None]

        entry_units = []
        for pu in proto:
            if len(units) >= len(CODES):
                break
            members = pu['members']
            uid = f"U{len(units) + 1:02d}"
            code = CODES[len(units)]
            # CA trace in the reference frame. Residues are numbered in order, the
            # second chain from 1001, so a viewer breaks the trace between chains.
            for k, m in enumerate(members):
                n = 0
                for key in sorted(ch[m]):
                    comp, ats = ch[m][key]
                    if 'CA' not in ats:
                        continue
                    n += 1
                    x, y, z = apply_rt(R, t, ats['CA'])
                    serial += 1
                    trace_lines.append(
                        f"ATOM  {serial:>5d}  CA  {(comp if comp in AA3 else 'GLY'):>3.3} "
                        f"{code:1}{k * 1000 + n:>4d}    {x:>8.3f}{y:>8.3f}{z:>8.3f}"
                        f"  1.00{min(99.99, rmsd):>6.2f}           C  ")
            units.append({
                'id': uid, 'pdb': pdb, 'chains': members,
                'chains_touching_antigen': [m for m in members if m in partners],
                'name': short_name(pu['descs']), 'descriptions': pu['descs'],
                'format': pu['format'], 'observed_lengths': pu['obs'], 'seqs': pu['seqs'],
                'resolution': meta.get('resolution'), 'rmsd': round(rmsd, 2),
                'epitope': sorted(pu['pairs']), 'n_pairs': pu['n_pairs'],
                'pairs': {str(r): n for r, n in sorted(pu['pairs'].items())},
                'near': {str(r): len(sx) for r, sx in sorted(pu['near'].items())},
                'chain_code': code,
            })
            entry_units.append(uid)

        ic_skipped = sum(1 for k in ch[c] if k[1] != '')
        entries.append({
            'pdb': pdb, 'status': 'ok', 'title': meta.get('title', ''),
            'resolution': meta.get('resolution'), 'method': meta.get('method'),
            'antigen_chain': c, 'offset': o, 'shift_to_reference': shift,
            'numbering': ('same as the reference' if shift == 0 else
                          'UniProt numbering' if o == 0 else
                          f'construct, offset {shift:+d} from the reference'),
            'superposition_rmsd': round(rmsd, 2), 'shared_ca': len(sh_pairs),
            'dropped_partners': dropped, 'dropped_copies': gone, 'units': entry_units,
            'antigen_residues_with_insertion_codes_skipped': ic_skipped,
        })
        fm = ', '.join(f"{u['format']}[{','.join(u['chains'])}]" for u in units[-len(entry_units):])
        print(f"  {pdb}  chain {c} shift {shift:+4d}  rmsd {rmsd:4.2f}  {fm}")

    ok = [e for e in entries if e.get('status') == 'ok']

    # ---- distinct binders ---------------------------------------------------
    print(f"\ncomparing {len(units)} units for sequence identity ...")
    parent = {u['id']: u['id'] for u in units}
    def f2(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    ident, similar = {}, []
    for i, ua in enumerate(units):
        for ub in units[i + 1:]:
            s = unit_identity(ua, ub)
            ident[(ua['id'], ub['id'])] = s
            if s >= args.identity:
                parent[f2(ua['id'])] = f2(ub['id'])
            elif s >= 0.80:
                similar.append({'a': ua['id'], 'b': ub['id'], 'identity': round(s, 3)})
    byid = {u['id']: u for u in units}
    groups = defaultdict(list)
    for u in units:
        groups[f2(u['id'])].append(u['id'])

    # the named binders: the cluster of the unit with the most contacts in that entry
    named = {}
    for nm, pdb in named_in.items():
        cand = [u for u in units if u['pdb'] == pdb]
        if not cand:
            print(f"  ! --name {nm}={pdb}: no unit found in that entry")
            continue
        named[nm] = max(cand, key=lambda u: u['n_pairs'])['id']

    clusters = []
    for k, (root, mem) in enumerate(sorted(groups.items(), key=lambda kv: min(kv[1])), 1):
        nm = next((n for n, uid in named.items() if f2(uid) == root), None)
        if nm:
            rep = named[nm]
        else:
            rep = min(mem, key=lambda m: (byid[m]['resolution'] or 99, byid[m]['rmsd']))
        clusters.append({'id': f'B{k:02d}', 'label': nm or byid[rep]['name'],
                         'display': (nm.capitalize() if nm else
                                     f"{byid[rep]['name'][:44]} ({byid[rep]['pdb']})"),
                         'named': bool(nm), 'rep': rep, 'format': byid[rep]['format'],
                         'members': sorted(mem), 'pdbs': sorted({byid[m]['pdb'] for m in mem})})
        for m in mem:
            byid[m]['cluster'] = f'B{k:02d}'; byid[m]['is_representative'] = (m == rep)

    # ---- lineage: does an antibody unit share a variable domain with a named binder? --
    for u in units:
        u['shares_chain_with'] = []
        if u['format'] not in ANTIBODY_LIKE:
            continue
        for nm, uid in named.items():
            ru = byid[uid]
            if ru['cluster'] == u['cluster']:
                continue
            best = max(nw_identity(sa[:VARIABLE_DOMAIN], sb[:VARIABLE_DOMAIN])
                       for sa in u['seqs'] for sb in ru['seqs'])
            if best >= args.identity:
                u['shares_chain_with'].append({'binder': nm, 'identity': round(best, 3)})

    # ---- footprints of the named binders, from all of their structures --------
    footprints = {}
    for cl in clusters:
        if not cl['named']:
            continue
        ref_u = byid[cl['rep']]
        allmem = [byid[m] for m in cl['members']]
        sims = {u['id']: (1.0 if u['id'] == ref_u['id'] else unit_identity(u, ref_u))
                for u in allmem}
        mem = [u for u in allmem if sims[u['id']] >= args.footprint_identity]
        variants = [{'unit': u['id'], 'pdb': u['pdb'], 'name': u['name'],
                     'identity': round(sims[u['id']], 3), 'epitope': u['epitope']}
                    for u in allmem if u not in mem]
        touched, modelled_n, pmax = defaultdict(int), defaultdict(int), defaultdict(int)
        for u in mem:
            for r in modelled_by_entry[u['pdb']]:
                modelled_n[r] += 1
            for r in u['epitope']:
                touched[r] += 1
                pmax[r] = max(pmax[r], u['pairs'][str(r)])
        union = sorted(touched)
        core = sorted(r for r in union if touched[r] / modelled_n[r] >= 0.5)
        footprints[cl['label']] = {
            'units': [u['id'] for u in mem], 'pdbs': sorted({u['pdb'] for u in mem}),
            'identity_floor': args.footprint_identity, 'variants': variants,
            'reference_unit': ref_u['id'], 'reference_entry': ref_u['pdb'],
            'reference_entry_only': ref_u['epitope'],
            'union': union, 'core': core,
            'support': {str(r): [touched[r], modelled_n[r]] for r in union},
            'pairs_max': {str(r): pmax[r] for r in union},
        }

    # ---- bins ---------------------------------------------------------------
    unions = {nm: set(fp['union']) for nm, fp in footprints.items()}
    for u in units:
        ep = set(u['epitope'])
        tally = defaultdict(int)
        for r in ep:
            tally[domain_of(r) or 'beyond'] += 1
        u['domain_counts'] = dict(tally)
        hit = {nm: len(ep & s) for nm, s in unions.items()}
        u['overlap'] = hit
        # a unit sits at a named binder's site when at least 3 residues AND a quarter of
        # its epitope lie in that footprint. Three residues out of sixteen is an edge.
        big = [nm for nm, n in hit.items() if n >= 3 and n >= 0.25 * len(ep)]
        if len(big) == 1:
            u['bin'] = f'{big[0]} site'
        elif len(big) > 1:
            u['bin'] = 'both named sites'
        elif domains and ep:
            top = max(tally, key=lambda d: tally[d])
            u['bin'] = f'domain {top}' if top != 'beyond' else 'beyond the reference'
        else:
            u['bin'] = 'other'

    # ---- counts, in the reference numbering ----------------------------------
    cu, cd, dens = defaultdict(int), defaultdict(int), defaultdict(int)
    for u in units:
        for r in u['epitope']:
            cu[r] += 1
    for cl in clusters:
        mem = [byid[m] for m in cl['members']]
        seen = {r for u in mem for r in u['epitope']}
        for r in seen:
            cd[r] += 1
        best = defaultdict(int)
        for u in mem:
            for r, n in u['near'].items():
                best[int(r)] = max(best[int(r)], n)
        for r, n in best.items():
            dens[r] += n
    top = max(ref_ca) if ref_ca else 0
    beyond = sorted(r for r in cu if r > top)

    for u in units:
        u.pop('seqs', None)                       # long, and not needed downstream
    skipped = defaultdict(int)
    for e in entries:
        if e.get('status') != 'ok':
            skipped[e['status']] += 1

    out = {
        'schema': 2,
        'generated': datetime.datetime.now().astimezone().isoformat(timespec='seconds'),
        'reference': os.path.basename(args.ref), 'reference_chain': args.ref_chain,
        'reference_offset': ref_off, 'reference_last_residue': top,
        'cutoff': args.cutoff, 'near_cutoff': args.near, 'identity_threshold': args.identity,
        'variable_domain_residues_compared': VARIABLE_DOMAIN,
        'domains': [{'name': n, 'first': a, 'last': b} for n, a, b in domains],
        'n_entries_searched': len(entries), 'n_entries_used': len(ok),
        'n_units': len(units), 'n_distinct': len(clusters),
        'skipped': dict(skipped),
        'sites': json.load(open(args.sites)) if args.sites else {},
        'footprints': footprints,
        'counts_units': {str(k): v for k, v in sorted(cu.items())},
        'counts_distinct': {str(k): v for k, v in sorted(cd.items())},
        'density_distinct': {str(k): v for k, v in sorted(dens.items())},
        'contacts_beyond_reference': beyond,
        'entries': entries, 'units': units, 'clusters': clusters,
        'similar_not_merged': sorted(similar, key=lambda s: -s['identity']),
    }
    with open(os.path.join(args.out, 'atlas.json'), 'w') as fh:
        json.dump(out, fh, indent=1)
        fh.write('\n')
    with open(os.path.join(args.out, 'binders.pdb'), 'w') as fh:
        fh.write('\n'.join(trace_lines) + '\nEND\n')

    # ---- report -------------------------------------------------------------
    print(f"\n{len(ok)} of {len(entries)} entries used -> {len(units)} units -> "
          f"{len(clusters)} distinct binders")
    for k in sorted(skipped, key=lambda k: -skipped[k]):
        print(f"  {skipped[k]:>3}  skipped: {k}")
    print("\ndistinct binders that appear more than once:")
    for cl in clusters:
        if len(cl['members']) > 1:
            print(f"  {cl['id']} {cl['label']:<28} {len(cl['members'])} units in "
                  f"{', '.join(cl['pdbs'])}")
    if similar:
        print("\nsimilar but not merged (identity 0.80 to threshold):")
        for s in sorted(similar, key=lambda s: -s['identity'])[:8]:
            a, b = byid[s['a']], byid[s['b']]
            print(f"  {s['identity']:.2f}  {a['pdb']} {a['name'][:24]:<24} | {b['pdb']} {b['name'][:24]}")
    for nm, fp in footprints.items():
        print(f"\nfootprint of {nm}: {len(fp['reference_entry_only'])} residues in "
              f"{fp['reference_entry']} alone, {len(fp['union'])} in {len(fp['pdbs'])} "
              f"entries ({', '.join(fp['pdbs'])}), {len(fp['core'])} in at least half")
        for v in fp['variants']:
            print(f"  not used (identity {v['identity']:.3f} to {nm}): {v['pdb']} {v['name'][:40]}")
        added = sorted(set(fp['union']) - set(fp['reference_entry_only']))
        if added:
            print(f"  added by the other entries: {runs_text(added)}")
    if beyond:
        print(f"\n{len(beyond)} residues beyond the reference ({top}) have contacts: {runs_text(beyond)}")
    print(f"\nwrote {args.out}/atlas.json and {args.out}/binders.pdb ({serial} CA atoms)")
    return 0


if __name__ == '__main__':
    sys.exit(main())
