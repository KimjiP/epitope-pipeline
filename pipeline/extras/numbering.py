import sys, math
from collections import defaultdict
AA3={'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
     'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
     'THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
cif, fas, CH = sys.argv[1], sys.argv[2], sys.argv[3]
lines=open(cif).read().splitlines()
cols,rows,i=[],[],0
while i<len(lines):
    if lines[i].startswith('loop_'):
        j=i+1;tags=[]
        while j<len(lines) and lines[j].startswith('_'): tags.append(lines[j].strip()); j+=1
        if tags and tags[0].startswith('_atom_site.'):
            cols=[t.split('.',1)[1] for t in tags];k=j
            while k<len(lines) and lines[k].strip() and not lines[k].startswith('#'):
                p=lines[k].split()
                if len(p)==len(cols): rows.append(p)
                k+=1
            i=k;continue
        i=j;continue
    i+=1
ix={c:n for n,c in enumerate(cols)}
seqres={}; nag=[]; ch_at=defaultdict(list)
for r in rows:
    ch=r[ix['auth_asym_id']]; comp=r[ix['label_comp_id']]
    xyz=(float(r[ix['Cartn_x']]),float(r[ix['Cartn_y']]),float(r[ix['Cartn_z']]))
    if r[ix['group_PDB']]=='ATOM':
        if ch==CH and comp in AA3: seqres[int(r[ix['auth_seq_id']])]=AA3[comp]
        if comp in AA3: ch_at[ch].append((int(r[ix['auth_seq_id']]),comp,r[ix['label_atom_id']],xyz))
    elif comp=='NAG':
        nag.append((ch,r[ix['auth_seq_id']],r[ix['label_atom_id']],xyz))

ref=''.join(l.strip() for l in open(fas) if not l.startswith('>'))
rr=sorted(seqres)
print(f"chain {CH}: first modelled residue = {rr[0]} ({seqres[rr[0]]}), last = {rr[-1]} ({seqres[rr[-1]]})")
probe=''.join(seqres[r] for r in rr[:30])
print(f"structure residues {rr[0]}-{rr[0]+29}: {probe}")
hit=ref.find(probe)
if hit>=0:
    off=(hit+1)-rr[0]
    print(f"\nFOUND in P04626 at UniProt position {hit+1}")
    print(f"OFFSET: UniProt number = structure number + {off}")
    for t in (557,558,560,569,571,572,573,592,593,603,310):
        if t in seqres: print(f"   structure {seqres[t]}{t}  =  UniProt {seqres[t]}{t+off}")
        else:           print(f"   structure {t} not modelled      UniProt {t+off} = {ref[t+off-1]}")
    print(f"\n   UniProt S310 (the famous mutation) = structure position {310-off} = {seqres.get(310-off,'?')}{310-off}")
else:
    print("probe not found - check chain")

print("\n=== observed NAG glycans: nearest Asn of chain", CH, "===")
asn=[(rs,at,p) for rs,comp,at,p in ch_at[CH] if comp=='ASN' and at=='ND2']
seen=set()
for ch,rs,at,p in nag:
    best=min(asn,key=lambda a: math.dist(a[2],p))
    d=math.dist(best[2],p)
    key=(ch,rs)
    if key not in seen and d<6:
        seen.add(key); print(f"  NAG {ch}{rs}: nearest ND2 = N{best[0]} at {d:.2f} A")
