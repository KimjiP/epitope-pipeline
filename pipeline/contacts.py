import sys, math
from collections import defaultdict

path = sys.argv[1]
CUT = float(sys.argv[2]) if len(sys.argv) > 2 else 4.5

lines = open(path).read().splitlines()
cols, atoms = [], []
i = 0
while i < len(lines):
    if lines[i].startswith('loop_'):
        j = i + 1
        tags = []
        while j < len(lines) and lines[j].startswith('_'):
            tags.append(lines[j].strip()); j += 1
        if tags and tags[0].startswith('_atom_site.'):
            cols = [t.split('.', 1)[1] for t in tags]
            k = j
            while k < len(lines) and lines[k].strip() and not lines[k].startswith('#'):
                p = lines[k].split()
                if len(p) == len(cols):
                    atoms.append(p)
                k += 1
            i = k; continue
        i = j; continue
    i += 1

ix = {c: n for n, c in enumerate(cols)}
def g(a, c): return a[ix[c]]

# report chains
chains = defaultdict(set)
for a in atoms:
    if g(a, 'group_PDB') != 'ATOM': continue
    chains[g(a, 'auth_asym_id')].add(int(g(a, 'auth_seq_id')))
print("CHAINS (auth id : n residues : first-last)")
for c in sorted(chains):
    r = sorted(chains[c])
    print(f"  {c} : {len(r):4d} : {r[0]}-{r[-1]}")

ANTIGEN = sys.argv[3] if len(sys.argv) > 3 else 'A'
PARTNER = list(sys.argv[4]) if len(sys.argv) > 4 else [c for c in chains if c != ANTIGEN]
print(f"\nantigen chain {ANTIGEN}  vs  partner chains {PARTNER}   cutoff {CUT} A")

ag, ab = [], []
for a in atoms:
    if g(a, 'group_PDB') != 'ATOM': continue
    el = g(a, 'type_symbol')
    if el == 'H': continue
    ch = g(a, 'auth_asym_id')
    xyz = (float(g(a, 'Cartn_x')), float(g(a, 'Cartn_y')), float(g(a, 'Cartn_z')))
    if ch == ANTIGEN:
        ag.append((int(g(a, 'auth_seq_id')), g(a, 'label_comp_id'), xyz))
    elif ch in PARTNER:
        ab.append((ch, int(g(a, 'auth_seq_id')), xyz))

# grid hash over antibody atoms
cell = CUT
grid = defaultdict(list)
for n, (ch, rs, (x, y, z)) in enumerate(ab):
    grid[(int(x//cell), int(y//cell), int(z//cell))].append(n)

contacts = defaultdict(int)
names = {}
nearest = {}
partners = defaultdict(set)
for rs, cn, (x, y, z) in ag:
    names[rs] = cn
    cx, cy, cz = int(x//cell), int(y//cell), int(z//cell)
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for dz in (-1, 0, 1):
                for n in grid.get((cx+dx, cy+dy, cz+dz), ()):
                    ch2, rs2, (x2, y2, z2) = ab[n]
                    d = math.dist((x, y, z), (x2, y2, z2))
                    if d <= CUT:
                        contacts[rs] += 1
                        partners[rs].add(f"{ch2}{rs2}")
                        if rs not in nearest or d < nearest[rs]:
                            nearest[rs] = d

AA3 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
       'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
       'THR':'T','TRP':'W','TYR':'Y','VAL':'V'}

print(f"\nEPITOPE: {len(contacts)} residues of chain {ANTIGEN} within {CUT} A of the partner\n")
print(f"{'resid':>6} {'aa':>4} {'1':>2} {'atom pairs':>11} {'closest A':>10}  touches")
for rs in sorted(contacts, key=lambda r: -contacts[r]):
    print(f"{rs:>6} {names[rs]:>4} {AA3.get(names[rs],'?'):>2} {contacts[rs]:>11} "
          f"{nearest[rs]:>10.2f}  {','.join(sorted(partners[rs]))}")

rr = sorted(contacts)
runs, start, prev = [], rr[0], rr[0]
for r in rr[1:]:
    if r - prev > 2:
        runs.append((start, prev)); start = r
    prev = r
runs.append((start, prev))
print("\nSEGMENTS (gaps of <=2 joined):")
for a2, b2 in runs:
    print(f"  {a2}-{b2}" if a2 != b2 else f"  {a2}")
print("\nhotspot shortlist (>=15 atom pairs):")
hs = [r for r in rr if contacts[r] >= 15]
print("  " + ",".join(f"{AA3.get(names[r],'X')}{r}" for r in hs))
print("  contig/hotspot style: " + ",".join(f"{ANTIGEN}{r}" for r in hs))
