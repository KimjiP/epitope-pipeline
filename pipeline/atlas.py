#!/usr/bin/env python3
"""Build an epitope atlas: every distinct binder of one antigen, in one frame.

Usage:
    python3 atlas.py REF.cif REF_CHAIN REFERENCE.fasta CACHE_DIR
        [--exclude FASTA] [--entities JSON] [--name NAME=PDB ...]
        [--cutoff 4.5] [--near 5.0] [--identity 0.95] [--min-cover 50] [--max-resolution 4.0]
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
4. Splits the binder chains into UNITS: one Fab (a heavy chain with its light
   chain), one nanobody, one scFv, or one non-antibody binder. Antibody chains
   are found from their sequence (the J-region motif that ends a variable
   domain), and two antibodies are never fused into one unit. A complex with
   two Fabs gives two units. Counting entries instead of units counts one
   antibody once for every entry it appears in.
5. Superposes the antigen onto the reference (Horn's quaternion method) and
   moves every unit into the reference frame. Each binder is placed with a fit
   on the antigen residues near ITS contacts (--local-radius), because one fit
   over the whole chain misplaces binders when the antigen moved as a whole.
6. Groups units into DISTINCT BINDERS. Antibodies: every chain type seen in
   both units must have a variable domain at least --identity the same (the
   lower of heavy and light counts), whatever the format. A unit that shows one
   chain only joins a complete antibody with the same CDR3. Other binders are
   compared over the whole chain. Entries worse than --max-resolution are
   skipped: side-chain contacts mean little there.
7. Counts, in the reference numbering, for every antigen residue:
      counts_units     how many units touch it
      counts_distinct  how many distinct binders touch it. A binder touches a
                       residue when it does so in at least half of its structures.
      counts_distinct_union  the same, counting any one structure
      density_distinct antibody residues within --near of it, summed over
                       distinct binders (the most any one structure shows)
8. With --name, builds a footprint for a named binder from ALL of its
   structures. One structure can miss a loop that another one models.

Everything is in the reference's numbering, so entries that use different
numbering conventions can be compared directly.

Stdlib only.
"""
import sys, os, re, json, math, glob, argparse, random, string, datetime
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


def overlap_identity(a, b, match=2, mismatch=-1, gap=-2):
    """Identity over the SHORTER of two sequences, aligned end to end with free end gaps.

    A copy that starts or stops a few residues sooner than another scores the same as the
    full one. nw_identity divides by the longer sequence and so punishes that; this one is
    for antibody variable domains, whose ends differ from file to file."""
    n, m = len(a), len(b)
    if n == 0 or m == 0:
        return 0.0
    H = [[0] * (m + 1) for _ in range(n + 1)]
    P = [[0] * (m + 1) for _ in range(n + 1)]       # 0 diagonal, 1 gap in b, 2 gap in a
    for i in range(1, n + 1):
        ai, Hi, Hp, Pi = a[i - 1], H[i], H[i - 1], P[i]
        for j in range(1, m + 1):
            s = Hp[j - 1] + (match if ai == b[j - 1] else mismatch)
            u = Hp[j] + gap
            l = Hi[j - 1] + gap
            if s >= u and s >= l:
                Hi[j] = s; Pi[j] = 0
            elif u >= l:
                Hi[j] = u; Pi[j] = 1
            else:
                Hi[j] = l; Pi[j] = 2
    bi, bj, best = n, 0, H[n][0]                    # free trailing gaps: best cell of the last row or column
    for j in range(m + 1):
        if H[n][j] > best:
            bi, bj, best = n, j, H[n][j]
    for i in range(n + 1):
        if H[i][m] > best:
            bi, bj, best = i, m, H[i][m]
    i, j, hits = bi, bj, 0
    while i > 0 and j > 0:
        p = P[i][j]
        if p == 0:
            if a[i - 1] == b[j - 1]:
                hits += 1
            i -= 1; j -= 1
        elif p == 1:
            i -= 1
        else:
            j -= 1
    return hits / min(n, m)


# ------------------------------------------------- antibody variable domains --

# A variable domain ends in its J region: WGQGTLVTVSS in a heavy-type domain (VH, VHH),
# FGGGTKLEIK in a light-type one. These motifs are conserved across species and germlines,
# so they find a domain without a numbering program. Only the conserved core is required: the
# last residues of the J region are often not modelled. A chain whose J region is missing
# altogether is not recognised; the name-based fallback in classify() then applies.
VH_J = re.compile(r'W[GA].G[TSAQ]..?[VL][TI]V[SF]?[SA]?')
VL_J = re.compile(r'F[GA].G[TS]...[VLIMF][KRL]?')


def v_domains(seq):
    """Variable domains in one chain, left to right: [(kind, vseq, cdr3)].

    kind is 'H' (heavy-type: VH or VHH) or 'L' (light-type). cdr3 is the stretch between the
    last Cys and the J region, which is enough to tell two clones apart. A chain that holds
    two or more domains is an scFv or a multivalent nanobody."""
    out, pos = [], 0
    while len(seq) - pos >= 95:
        mh, ml = VH_J.search(seq, pos + 85), VL_J.search(seq, pos + 85)
        cands = [(m.end(), k, m) for k, m in (('H', mh), ('L', ml)) if m]
        if not cands:
            break
        end, kind, m = min(cands, key=lambda c: c[0])
        start = pos
        if out:                                     # a later domain: skip the Gly/Ser linker before it
            while start < m.start() and seq[start] in 'GS':
                start += 1
        v = seq[start:end]
        j = m.start() - start
        c = v.rfind('C', 0, j)
        cdr3 = v[c + 1:j] if c >= 0 and 4 <= j - c <= 40 else ''
        out.append((kind, v, cdr3))
        pos = end
    return out


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


SCAFFOLD_WORDS = (('ankyrin', 'DARPin'), ('darpin', 'DARPin'), ('repebody', 'Repebody'),
                  ('affibody', 'Affibody'), ('immunoglobulin g-binding protein a', 'Affibody'),
                  ('adnectin', 'Adnectin'), ('monobody', 'Monobody'), ('lipocalin', 'Anticalin'),
                  ('anticalin', 'Anticalin'), ('minibinder', 'Designed binder'),
                  ('de novo', 'Designed binder'))


def classify(descs, title, obs_lens, v=None, one_chain_both=False):
    """(kind, format) of one binding unit.

    kind  'antibody'  a variable domain was found in the sequence (or, failing that, the name says so)
          'scaffold'  an engineered non-antibody binder: DARPin, repebody, affibody, Fcab, ...
          'other'     anything else. Often the antigen's own ligand or receptor.
    The sequence decides first. Names in the PDB are free text and often say only 'heavy chain'."""
    n = len(obs_lens)
    d = ' '.join(descs).lower()
    t = (title or '').lower()
    if v:
        if one_chain_both:
            return 'antibody', 'scFv'
        if v.get('H') and v.get('L'):
            return 'antibody', ('Fab' if max(obs_lens) > 150 else 'Fv')
        if v.get('H'):
            multi = len(v['H']) > 1
            return 'antibody', ('VHH' if multi or max(obs_lens) <= 160 else 'antibody fragment')
        return 'antibody', 'antibody fragment'
    for word, fmt in SCAFFOLD_WORDS:
        if word in d:
            return 'scaffold', fmt
    if 'ig gamma' in d and n == 1:
        return 'scaffold', 'Fcab (engineered Fc)'
    if 'nanobod' in d or 'vhh' in d:
        return 'antibody', 'VHH'
    if 'scfv' in d or 'single-chain' in d or 'scfv' in t:
        return 'antibody', 'scFv'
    if any(w in d for w in ('fab', 'antibody', 'immunoglobulin', 'igg')) or ('heavy' in d and 'light' in d):
        return 'antibody', ('Fab' if n == 2 else 'antibody fragment')
    return 'other', 'other protein'


def comparable_seq(seq, fmt):
    """The part of a chain that carries the binding loops. An Fcab (engineered Fc) is compared
    over its CH3 domain, where the loops were engineered: the CH2 domain is the same in every
    Fcab, and over the whole Fc two different Fcabs look 97% alike. Other binders: the whole chain."""
    return seq[-110:] if fmt.startswith('Fcab') else seq


def unit_identity(ua, ub):
    """How alike two binding units are, 0 to 1.

    Antibodies: every chain type (heavy, light) that BOTH units show must match. The score is the
    LOWER of the variable-domain identities, so a shared light chain cannot hide a different heavy
    chain. The format, the number of chains and the constant domains play no part: a Fab and an
    scFv of one clone are one binder.
    Other binders: the mean identity of the chains; 0 when the format or the chain count differ."""
    if ua['kind'] != ub['kind']:
        return 0.0
    va, vb = ua.get('v'), ub.get('v')
    if va and vb:
        comps = [max(overlap_identity(x[0], y[0]) for x in va[t] for y in vb[t])
                 for t in ('H', 'L') if va.get(t) and vb.get(t)]
        return min(comps) if comps else 0.0
    if va or vb:
        return 0.0
    if ua['format'] != ub['format'] or len(ua['seqs']) != len(ub['seqs']):
        return 0.0
    A = [comparable_seq(s, ua['format']) for s in ua['seqs']]
    B = [comparable_seq(s, ub['format']) for s in ub['seqs']]
    for a, b in zip(sorted(A, key=len), sorted(B, key=len)):
        if min(len(a), len(b)) < 0.6 * max(len(a), len(b)):
            return 0.0
    # identity over the shorter chain: the same protein modelled a few residues longer or shorter
    # (a different construct, a different disordered tail) is still the same binder
    if len(A) == 1:
        return overlap_identity(A[0], B[0])
    if len(A) == 2:
        straight = (overlap_identity(A[0], B[0]) + overlap_identity(A[1], B[1])) / 2
        crossed = (overlap_identity(A[0], B[1]) + overlap_identity(A[1], B[0])) / 2
        return max(straight, crossed)
    return sum(overlap_identity(a, b) for a, b in zip(A, B)) / len(A)


def cluster_units(units, thr):
    """Group binding units into distinct binders. Returns (groups, near_misses).

    1. Complete antibodies (a heavy and a light variable domain seen) and all non-antibodies are
       clustered with COMPLETE linkage: every pair in a group must reach the threshold. Single
       linkage let one in-between unit chain two different clones into one binder.
    2. A unit that shows only one chain (a VH with no modelled VL, a nanobody) joins the closest
       complete antibody when its chain matches and its CDR3 is identical. It never joins two
       groups together.
    3. Single-chain units that match no complete antibody (nanobodies) are clustered among
       themselves."""
    byid = {u['id']: u for u in units}
    memo = {}

    def score(a, b):
        k = (a, b) if a < b else (b, a)
        if k not in memo:
            memo[k] = unit_identity(byid[a], byid[b])
        return memo[k]

    def full(u):
        return bool(u.get('v') and u['v'].get('H') and u['v'].get('L'))

    def partial(u):
        return bool(u.get('v')) and not full(u)

    def complete_link(ids):
        parent = {i: i for i in ids}
        members = {i: [i] for i in ids}

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]; x = parent[x]
            return x
        pairs = sorted(((score(a, b), a, b) for n, a in enumerate(ids) for b in ids[n + 1:]),
                       key=lambda t: (-t[0], t[1], t[2]))
        for s, a, b in pairs:
            if s < thr:
                break
            ra, rb = find(a), find(b)
            if ra != rb and all(score(x, y) >= thr for x in members[ra] for y in members[rb]):
                parent[ra] = rb
                members[rb] += members.pop(ra)
        return list(members.values())

    groups = complete_link([u['id'] for u in units if not partial(u)])
    loose = []
    for u in (x for x in units if partial(x)):
        t = 'H' if u['v'].get('H') else 'L'
        best = None
        for gi, g in enumerate(groups):
            for m in g:
                vm = byid[m].get('v')
                for x in u['v'][t]:
                    for y in (vm or {}).get(t, []):
                        s = overlap_identity(x[0], y[0])
                        if s >= thr and (not x[1] or not y[1] or x[1] == y[1]):
                            key = (s, len(g), -gi)
                            if best is None or key > best[0]:
                                best = (key, gi)
        if best:
            groups[best[1]].append(u['id'])
        else:
            loose.append(u['id'])
    groups += complete_link(loose)
    similar = [{'a': a, 'b': b, 'identity': round(s, 3)} for (a, b), s in memo.items() if 0.80 <= s < thr]
    return groups, similar


def local_fit(ent_ca, ref_ca, epitope, radius, fallback):
    """Fit the antigen onto the reference on the rigid core around the binder's contacts.

    One fit over the whole chain puts a binder in the wrong place when the antigen as a whole has
    moved (EGFR, tethered against extended: up to 25 A). So the fit uses the antigen residues
    within `radius` of the binder's contacts, and then drops, cycle by cycle, the residues that
    sit more than 2 A (or twice the rmsd) from their partner, as structure-alignment programs do.
    What is left is the domain the binder sits on. Returns (R, t, rmsd, n_residues, scope);
    `fallback` is the whole-chain fit, used when too few residues are near the contacts."""
    core = [ent_ca[r] for r in epitope if r in ent_ca]
    if not core:
        return tuple(fallback) + (0, 'global')
    keep = [r for r, p in ent_ca.items()
            if r in ref_ca and any(math.dist(p, q) <= radius for q in core)]
    if len(keep) < 12:
        return tuple(fallback) + (len(keep), 'global')
    R, t, rmsd = superpose([ent_ca[r] for r in keep], [ref_ca[r] for r in keep])
    for _ in range(5):
        dev = {r: math.dist(apply_rt(R, t, ent_ca[r]), ref_ca[r]) for r in keep}
        kept = [r for r in keep if dev[r] <= max(2.0, 2.0 * rmsd)]
        if len(kept) < 12 or len(kept) == len(keep):
            break
        keep = kept
        R, t, rmsd = superpose([ent_ca[r] for r in keep], [ref_ca[r] for r in keep])
    return R, t, rmsd, len(keep), 'local'



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
                              'released': ((e.get('rcsb_accession_info') or {})
                                           .get('initial_release_date') or '')[:10],
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
    ap.add_argument('--max-resolution', type=float, default=0.0, metavar='A',
                    help='skip entries solved at a resolution worse than this, in A (0 = keep all). '
                         'Side-chain contacts at 4.5 A mean little in a 7 A map.')
    ap.add_argument('--local-radius', type=float, default=15.0, metavar='A',
                    help='each binder is placed with a fit on the antigen residues this close to its contacts')
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
        res_here = meta.get('resolution')
        if args.max_resolution and res_here and res_here > args.max_resolution:
            entries.append({'pdb': pdb, 'resolution': res_here,
                            'status': f'resolution worse than {args.max_resolution:g} A'})
            continue

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

        # Units. A Fab is a heavy chain with its light chain; an scFv or a nanobody is one chain.
        # Chains are typed from their sequence, and each heavy chain is paired with the light chain
        # it shares the most atom pairs with. Two antibodies are NEVER fused. (Fusing every pair of
        # touching chains put a Fab together with its crystal-packing neighbour, a 3 to 6 chain
        # "unit" that no other entry could match, and one drug was counted three times.)
        # Chains with no variable domain are linked by contact, as before.
        pts = {n: [p for _, p in cand[n]] for n in cand}
        cseq = {n: seq_of(ch[n]) for n in cand}
        cdom = {n: v_domains(cseq[n]) for n in cand}
        ctype = {n: (None if not d else ('S' if len(d) > 1 else d[0][0])) for n, d in cdom.items()}
        hs = [n for n in sorted(cand) if ctype[n] == 'H']
        ls = [n for n in sorted(cand) if ctype[n] == 'L']
        hl = sorted(((chain_pairs(pts[h], pts[l], args.cutoff, 10 ** 9), h, l) for h in hs for l in ls),
                    key=lambda t: (-t[0], t[1], t[2]))
        used, glist = set(), []
        for k, h, l in hl:
            if k >= args.link_pairs and h not in used and l not in used:
                used |= {h, l}
                glist.append([h, l])
        glist += [[n] for n in sorted(cand) if ctype[n] is not None and n not in used]
        plain = [n for n in sorted(cand) if ctype[n] is None]
        parent = {n: n for n in plain}

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]; x = parent[x]
            return x
        for i, a in enumerate(plain):
            for b in plain[i + 1:]:
                if a in partners or b in partners:
                    if chain_pairs(pts[a], pts[b], args.cutoff, args.link_pairs) >= args.link_pairs:
                        parent[find(a)] = find(b)
        by_root = defaultdict(list)
        for n in plain:
            by_root[find(n)].append(n)
        glist += list(by_root.values())
        groups = [g for g in glist if any(m in partners for m in g)]

        ent_ca = {r + shift: ats['CA'] for r, (comp, ats) in agres.items() if 'CA' in ats}
        sh_pairs = [(ent_ca[r], ref_ca[r]) for r in ent_ca if r in ref_ca]
        if len(sh_pairs) < 30:
            entries.append({'pdb': pdb, 'status': 'too few shared CA atoms to superpose'})
            continue
        R, t, rmsd = superpose([p[0] for p in sh_pairs], [p[1] for p in sh_pairs])
        modelled_by_entry[pdb] = set(modelled)
        # other copies of the antigen in this entry: contacts to them are not part of the epitope
        # (they may be crystal packing), but they are counted so that a quaternary epitope shows up
        other_atoms = [(0, p) for c3 in ag if c3 != c
                       for _, ats in intkeyed(ch[c3]).values() for p in ats.values()]

        proto = []
        for members in sorted(groups, key=lambda g: min(g)):
            members = sorted(members, key=lambda x: -len(ch[x]))
            uatoms = [a for m in members for a in cand[m]]
            pairs, near = interface(ag_atoms, uatoms, args.cutoff, args.near)
            descs = [meta.get('desc', {}).get(m, '') for m in members]
            obs = [len(ch[m]) for m in members]
            v = defaultdict(list)
            for m in members:
                for k, vs, c3 in cdom[m]:
                    v[k].append((vs, c3))
            both = any(len({k for k, _, _ in cdom[m]}) == 2 for m in members)
            kind, fmt = classify(descs, meta.get('title'), obs, dict(v) or None, both)
            other_pairs = (sum(interface(other_atoms, uatoms, args.cutoff, args.cutoff)[0].values())
                           if other_atoms else 0)
            proto.append({'members': members, 'pairs': pairs, 'near': near, 'descs': descs,
                          'obs': obs, 'kind': kind, 'format': fmt, 'v': dict(v) or None,
                          'seqs': [cseq[m] for m in members], 'n_pairs': sum(pairs.values()),
                          'other_pairs': other_pairs})

        # Identical copies in one entry: only the copy that really binds is a binder. Index by
        # position and re-read the list each time. Testing the loop variable instead let an entry
        # with three or more identical copies drop the same one twice, which raised ValueError on an
        # EGFR entry carrying four copies of nanobody 7D12.
        gone = []
        for i in range(len(proto)):
            for j in range(i + 1, len(proto)):
                a, b = proto[i], proto[j]
                if a is None or b is None or unit_identity(a, b) < args.identity:
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
            # Each binder is placed with a fit on the antigen residues near ITS contacts, not on the
            # whole chain. A CA trace in the reference frame. Residues are numbered in order, the
            # second chain from 1001, so a viewer breaks the trace between chains.
            Ru, tu, rmsd_u, n_fit, scope = local_fit(ent_ca, ref_ca, sorted(pu['pairs']),
                                                      args.local_radius, (R, t, rmsd))
            for k, m in enumerate(members):
                n = 0
                for key in sorted(ch[m]):
                    comp, ats = ch[m][key]
                    if 'CA' not in ats:
                        continue
                    n += 1
                    x, y, z = apply_rt(Ru, tu, ats['CA'])
                    serial += 1
                    trace_lines.append(
                        f"ATOM  {serial:>5d}  CA  {(comp if comp in AA3 else 'GLY'):>3.3} "
                        f"{code:1}{k * 1000 + n:>4d}    {x:>8.3f}{y:>8.3f}{z:>8.3f}"
                        f"  1.00{min(99.99, rmsd_u):>6.2f}           C  ")
            units.append({
                'id': uid, 'pdb': pdb, 'chains': members,
                'chains_touching_antigen': [m for m in members if m in partners],
                'name': short_name(pu['descs']), 'descriptions': pu['descs'],
                'kind': pu['kind'], 'format': pu['format'], 'observed_lengths': pu['obs'],
                'v': pu['v'], 'seqs': pu['seqs'],
                'resolution': meta.get('resolution'), 'rmsd': round(rmsd_u, 2),
                'rmsd_global': round(rmsd, 2), 'fit': scope, 'fit_residues': n_fit,
                'epitope': sorted(pu['pairs']), 'n_pairs': pu['n_pairs'],
                'other_copy_pairs': pu['other_pairs'],
                'other_copy_share': round(pu['other_pairs'] / max(1, pu['other_pairs'] + pu['n_pairs']), 3),
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
    glist, similar = cluster_units(units, args.identity)
    byid = {u['id']: u for u in units}
    groups = {min(g): sorted(g) for g in glist}
    root_of = {m: root for root, g in groups.items() for m in g}

    # the named binders: the cluster of the unit with the most contacts in that entry
    named = {}
    for nm, pdb in named_in.items():
        cand = [u for u in units if u['pdb'] == pdb]
        if not cand:
            print(f"  ! --name {nm}={pdb}: no unit found in that entry")
            continue
        named[nm] = max(cand, key=lambda u: u['n_pairs'])['id']

    clusters = []
    for k, (root, mem) in enumerate(sorted(groups.items()), 1):
        nm = next((n for n, uid in named.items() if root_of[uid] == root), None)
        if nm:
            rep = named[nm]
        else:
            rep = min(mem, key=lambda m: (byid[m]['resolution'] or 99, byid[m]['rmsd']))
        clusters.append({'id': f'B{k:02d}', 'label': nm or byid[rep]['name'],
                         'display': (nm.capitalize() if nm else
                                     f"{byid[rep]['name'][:44]} ({byid[rep]['pdb']})"),
                         'named': bool(nm), 'rep': rep, 'kind': byid[rep]['kind'],
                         'format': byid[rep]['format'],
                         'members': sorted(mem), 'pdbs': sorted({byid[m]['pdb'] for m in mem})})
        for m in mem:
            byid[m]['cluster'] = f'B{k:02d}'; byid[m]['is_representative'] = (m == rep)

    # ---- lineage: does an antibody unit share a variable domain with a named binder? --
    for u in units:
        u['shares_chain_with'] = []
        if u['kind'] != 'antibody' or not u.get('v'):
            continue
        for nm, uid in named.items():
            ru = byid[uid]
            if ru['cluster'] == u['cluster'] or not ru.get('v'):
                continue
            best = max((overlap_identity(x[0], y[0]) for t in ('H', 'L')
                        for x in u['v'].get(t, []) for y in ru['v'].get(t, [])), default=0.0)
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
    # A binder TOUCHES a residue when it does so in at least half of its structures (in all of them
    # when it was solved once or twice). The union over every structure gave a popular drug a larger
    # footprint than any single structure shows: trastuzumab, 17 to 20 residues per structure, 28 in
    # the union. counts_distinct_union keeps the inclusive version.
    cu, cd, cdu, dens = defaultdict(int), defaultdict(int), defaultdict(int), defaultdict(int)
    for u in units:
        for r in u['epitope']:
            cu[r] += 1
    for cl in clusters:
        mem = [byid[m] for m in cl['members']]
        tally = defaultdict(int)
        for u in mem:
            for r in u['epitope']:
                tally[r] += 1
        need = (len(mem) + 1) // 2
        cl['epitope_core'] = sorted(r for r, k in tally.items() if k >= need)
        cl['epitope_union'] = sorted(tally)
        for r in cl['epitope_core']:
            cd[r] += 1
        for r in cl['epitope_union']:
            cdu[r] += 1
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
        'schema': 3,
        'generated': datetime.datetime.now().astimezone().isoformat(timespec='seconds'),
        'reference': os.path.basename(args.ref), 'reference_chain': args.ref_chain,
        'reference_offset': ref_off, 'reference_last_residue': top,
        'cutoff': args.cutoff, 'near_cutoff': args.near, 'identity_threshold': args.identity,
        'max_resolution': args.max_resolution or None, 'local_radius': args.local_radius,
        'domains': [{'name': n, 'first': a, 'last': b} for n, a, b in domains],
        'n_entries_searched': len(entries), 'n_entries_used': len(ok),
        'n_units': len(units), 'n_distinct': len(clusters),
        'n_antibodies': sum(1 for c in clusters if c['kind'] == 'antibody'),
        'skipped': dict(skipped),
        'sites': json.load(open(args.sites)) if args.sites else {},
        'footprints': footprints,
        'counts_units': {str(k): v for k, v in sorted(cu.items())},
        'counts_distinct': {str(k): v for k, v in sorted(cd.items())},
        'counts_distinct_union': {str(k): v for k, v in sorted(cdu.items())},
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
