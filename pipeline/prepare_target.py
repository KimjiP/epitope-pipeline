#!/usr/bin/env python3
"""Prepare a design-model target: a clean structure file plus a hotspot list.

Usage
-----
    python3 prepare_target.py STRUCTURE CHAIN --window LO-HI \
        --hotspots E558,D560,P571 [--renumber] [--uniprot REF.fasta] \
        [--outdir DIR] [--name NAME]

Example
-------
    python3 prepare_target.py 1N8Z.cif C --window 490-607 \
        --hotspots E558,D560,P571,P572,F573,I591,W592,P603 \
        --uniprot P04626.fasta --name her2_domainIV

What it writes
--------------
    <name>_target.pdb       the chain, the window, nothing else
    <name>_hotspots.txt     the model-ready residue list
    <name>_numbering.json   every mapping, gap and digest

Why the hotspot spec carries an amino-acid letter
-------------------------------------------------
In the model's own notation, `C558` means CHAIN C, residue 558. In ordinary
protein notation, `C558` means CYSTEINE 558. The two collide.

So this script takes only the amino-acid form on input (`E558` = glutamate at
558) and VERIFIES it against the structure. A wrong letter stops the run. The
script then emits the chain-prefixed form that the model wants.

Renumbering
-----------
Design pipelines often renumber a structure continuously. That changes every
residue number, so a hotspot list written against the original file becomes
wrong and stays silent about it. `--renumber` does the renumbering here, remaps
the hotspots itself, and records both numbers side by side.

Stdlib only. Exits non-zero if any check fails.
"""
import sys, os, json, argparse, hashlib, datetime
from collections import defaultdict

AA3 = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q',
       'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
       'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
       'TYR': 'Y', 'VAL': 'V', 'MSE': 'M', 'SEC': 'U', 'PYL': 'O'}
BACKBONE = ('N', 'CA', 'C', 'O')


# ---------------------------------------------------------------- reading ----

def read_structure(path):
    """Return a list of atom dicts from an mmCIF or PDB file. ATOM records only."""
    text = open(path).read()
    if '_atom_site.' in text:
        return _read_cif(text)
    return _read_pdb(text)


def _read_cif(text):
    lines = text.splitlines()
    cols, rows, i = [], [], 0
    while i < len(lines):
        if lines[i].startswith('loop_'):
            j, tags = i + 1, []
            while j < len(lines) and lines[j].startswith('_'):
                tags.append(lines[j].strip())
                j += 1
            if tags and tags[0].startswith('_atom_site.'):
                cols = [t.split('.', 1)[1] for t in tags]
                k = j
                while k < len(lines) and lines[k].strip() and not lines[k].startswith('#'):
                    p = lines[k].split()
                    if len(p) == len(cols):
                        rows.append(p)
                    k += 1
                i = k
                continue
            i = j
            continue
        i += 1
    if not cols:
        sys.exit('no _atom_site loop found in the mmCIF file')
    ix = {c: n for n, c in enumerate(cols)}

    def get(r, *names):
        for n in names:
            if n in ix:
                v = r[ix[n]]
                if v not in ('.', '?'):
                    return v
        return ''

    out = []
    for r in rows:
        if r[ix['group_PDB']] != 'ATOM':
            continue
        alt = get(r, 'label_alt_id')
        if alt not in ('', 'A'):
            continue
        out.append(dict(
            name=get(r, 'label_atom_id', 'auth_atom_id'),
            resname=get(r, 'label_comp_id', 'auth_comp_id'),
            chain=get(r, 'auth_asym_id', 'label_asym_id'),
            resseq=int(get(r, 'auth_seq_id', 'label_seq_id')),
            icode=get(r, 'pdbx_PDB_ins_code'),
            x=float(r[ix['Cartn_x']]), y=float(r[ix['Cartn_y']]), z=float(r[ix['Cartn_z']]),
            occ=float(get(r, 'occupancy') or 1.0),
            bfac=float(get(r, 'B_iso_or_equiv') or 0.0),
            element=get(r, 'type_symbol')))
    return out


def _read_pdb(text):
    out = []
    for line in text.splitlines():
        if not line.startswith('ATOM'):
            continue
        alt = line[16].strip()
        if alt not in ('', 'A'):
            continue
        el = line[76:78].strip()
        name = line[12:16].strip()
        out.append(dict(
            name=name, resname=line[17:20].strip(), chain=line[21].strip() or ' ',
            resseq=int(line[22:26]), icode=line[26].strip(),
            x=float(line[30:38]), y=float(line[38:46]), z=float(line[46:54]),
            occ=float(line[54:60] or 1.0), bfac=float(line[60:66] or 0.0),
            element=el or (name[0] if name else 'C')))
    if not out:
        sys.exit('no ATOM records found in the PDB file')
    return out


# ---------------------------------------------------------------- writing ----

def pdb_line(serial, a):
    """One PDB ATOM record. Columns per the PDB 3.3 specification."""
    name, el = a['name'], (a['element'] or 'C')
    if len(name) >= 4 or len(el) == 2:
        nm = f'{name:<4.4}'
    else:
        nm = f' {name:<3.3}'
    return (f"ATOM  {serial:>5d} {nm}{'':1}{a['resname']:>3.3} {a['chain'][:1]:1}"
            f"{a['resseq']:>4d}{(a['icode'] or ' '):1}   "
            f"{a['x']:>8.3f}{a['y']:>8.3f}{a['z']:>8.3f}"
            f"{a['occ']:>6.2f}{a['bfac']:>6.2f}          {el:>2.2}  ")


def write_pdb(path, atoms, header_lines):
    with open(path, 'w') as fh:
        for h in header_lines:
            fh.write(f'REMARK 999 {h}\n'[:81].rstrip() + '\n')
        serial = 0
        for a in atoms:
            serial += 1
            fh.write(pdb_line(serial, a) + '\n')
        if atoms:
            last = atoms[-1]
            serial += 1
            fh.write(f"TER   {serial:>5d}      {last['resname']:>3.3} "
                     f"{last['chain'][:1]:1}{last['resseq']:>4d}\n")
        fh.write('END\n')


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for blk in iter(lambda: fh.read(65536), b''):
            h.update(blk)
    return h.hexdigest()


# ------------------------------------------------------------------ logic ----

def residue_table(atoms):
    """{resseq: (resname, one_letter)} and the ordered residue list."""
    seen = {}
    order = []
    for a in atoms:
        if a['resseq'] not in seen:
            seen[a['resseq']] = a['resname']
            order.append(a['resseq'])
    return {r: (seen[r], AA3.get(seen[r], 'X')) for r in order}, sorted(order)


def find_gaps(resids):
    return [(a + 1, b - 1) for a, b in zip(resids, resids[1:]) if b - a > 1]


def parse_hotspots(spec):
    """'E558,D560,571' -> [(expected_letter_or_None, number), ...]"""
    out = []
    for tok in spec.replace(' ', '').split(','):
        if not tok:
            continue
        letters = ''.join(c for c in tok if c.isalpha())
        digits = ''.join(c for c in tok if c.isdigit())
        if not digits:
            sys.exit(f'hotspot "{tok}" has no residue number')
        if len(letters) > 1:
            sys.exit(f'hotspot "{tok}": use one amino-acid letter, e.g. E558')
        out.append((letters.upper() or None, int(digits)))
    if not out:
        sys.exit('--hotspots is empty')
    return out


def uniprot_offset(ref_path, table, resids):
    """Align the first 30 modelled residues to a reference FASTA."""
    seq = ''.join(l.strip() for l in open(ref_path) if not l.startswith('>'))
    probe_ids = [r for r in resids if r == resids[0] or r - resids[0] < 40][:30]
    probe = ''.join(table[r][1] for r in probe_ids)
    if len(probe) < 10:
        return None, None, 'probe too short to align'
    hit = seq.find(probe)
    if hit < 0:
        return None, None, 'probe not found in the reference sequence'
    if seq.find(probe, hit + 1) != -1:
        return None, None, 'probe is not unique in the reference sequence'
    return (hit + 1) - probe_ids[0], len(seq), None


# ------------------------------------------------------------------- main ----

def main():
    ap = argparse.ArgumentParser(
        description='Prepare a design-model target file and a verified hotspot list.',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('structure', help='input .cif or .pdb')
    ap.add_argument('chain', help='antigen chain, as the input file labels it')
    ap.add_argument('--window', required=True, metavar='LO-HI',
                    help='residue range to keep, e.g. 490-607')
    ap.add_argument('--hotspots', required=True,
                    help='comma list in amino-acid form, e.g. E558,D560,P571')
    ap.add_argument('--renumber', action='store_true',
                    help='renumber the kept residues continuously from 1 and remap the hotspots')
    ap.add_argument('--uniprot', metavar='FASTA',
                    help='reference sequence, to record the numbering offset')
    ap.add_argument('--outdir', default='target', help='output directory (default: target)')
    ap.add_argument('--name', default='target', help='file name stem (default: target)')
    args = ap.parse_args()

    try:
        lo, hi = (int(v) for v in args.window.split('-'))
    except ValueError:
        sys.exit('--window must look like 490-607')
    if lo > hi:
        sys.exit('--window: LO is greater than HI')

    problems = []
    notes = []

    atoms = read_structure(args.structure)
    chains = sorted({a['chain'] for a in atoms})
    if args.chain not in chains:
        sys.exit(f'chain {args.chain} not in the file. Chains present: {chains}')

    full = [a for a in atoms if a['chain'] == args.chain]
    full_table, full_ids = residue_table(full)

    kept = [a for a in full if lo <= a['resseq'] <= hi and not a['icode']]
    if not kept:
        sys.exit(f'no residues of chain {args.chain} in {lo}-{hi}')
    dropped_icode = {a['resseq'] for a in full if a['icode'] and lo <= a['resseq'] <= hi}
    if dropped_icode:
        notes.append(f'dropped insertion-coded residues: {sorted(dropped_icode)}')

    table, resids = residue_table(kept)
    gaps = find_gaps(resids)

    # --- hotspots, checked against the structure --------------------------
    hotspots = parse_hotspots(args.hotspots)
    checked = []
    for want, num in hotspots:
        if num not in table:
            if num in full_table:
                problems.append(f'hotspot {num} is outside the window {lo}-{hi} '
                                f'(it exists in the chain as '
                                f'{full_table[num][1]}{num})')
            else:
                problems.append(f'hotspot {num} is not modelled in chain {args.chain}')
            continue
        got = table[num][1]
        if want and want != got:
            problems.append(f'hotspot {want}{num} is wrong: residue {num} is '
                            f'{got} ({table[num][0]}), not {want}')
            continue
        checked.append((got, num))

    # --- optional renumbering --------------------------------------------
    old_to_new = {r: i for i, r in enumerate(resids, 1)}
    if args.renumber:
        for a in kept:
            a['resseq'] = old_to_new[a['resseq']]
        notes.append('renumbered continuously from 1; the hotspot list below is remapped')
    emit = [(aa, old_to_new[n] if args.renumber else n) for aa, n in checked]

    # --- reference-sequence offset ---------------------------------------
    offset = seqlen = off_err = None
    if args.uniprot:
        offset, seqlen, off_err = uniprot_offset(args.uniprot, full_table, full_ids)
        if off_err:
            notes.append(f'reference alignment: {off_err}')

    os.makedirs(args.outdir, exist_ok=True)
    stem = os.path.join(args.outdir, args.name)
    pdb_path, hot_path, json_path = f'{stem}_target.pdb', f'{stem}_hotspots.txt', f'{stem}_numbering.json'

    if problems:
        print('FAILED - nothing written\n')
        for p in problems:
            print(f'  ! {p}')
        print('\nFix the hotspot list or the window, then run again.')
        return 2

    hot_model = ','.join(f'{args.chain}{n}' for _, n in emit)
    header = [
        f'prepared by prepare_target.py on {datetime.date.today().isoformat()}',
        f'source {os.path.basename(args.structure)} chain {args.chain} residues {lo}-{hi}',
        f'hotspots {hot_model}',
        'renumbered continuously from 1' if args.renumber else 'original numbering kept',
    ]
    write_pdb(pdb_path, kept, header)

    with open(hot_path, 'w') as fh:
        fh.write(f'# hotspot residues on chain {args.chain} of '
                 f'{os.path.basename(pdb_path)}\n')
        fh.write('# the leading letter is the CHAIN, not the amino acid\n')
        fh.write(f'{hot_model}\n')
        fh.write('\n# the same residues, as amino acid + number in this file\n')
        fh.write('# ' + ','.join(f'{aa}{n}' for aa, n in emit) + '\n')
        if offset is not None:
            fh.write('# the same residues in reference numbering '
                     f'(this file + {offset})\n')
            fh.write('# ' + ','.join(f'{a}{n + offset}' for a, n in checked) + '\n')

    record = {
        'generated': datetime.datetime.now().astimezone().isoformat(timespec='seconds'),
        'tool': 'prepare_target.py',
        'source': {'file': os.path.abspath(args.structure),
                   'sha256': sha256(args.structure),
                   'chains_present': chains,
                   'chain_used': args.chain},
        'window': {'requested': [lo, hi],
                   'kept_first_last': [resids[0], resids[-1]],
                   'residues_kept': len(resids),
                   'atoms_kept': len(kept)},
        'gaps_in_window': [list(g) for g in gaps],
        'renumbered': bool(args.renumber),
        'numbering': {
            'file_to_original': ({str(v): k for k, v in old_to_new.items()}
                                 if args.renumber else 'identity'),
            'original_to_reference_offset': offset,
            'reference_file': (os.path.abspath(args.uniprot) if args.uniprot else None),
            'reference_length': seqlen,
        },
        'hotspots': [{'amino_acid': aa,
                      'original': orig,
                      'in_file': old_to_new[orig] if args.renumber else orig,
                      'reference': (orig + offset) if offset is not None else None,
                      'model_syntax': f'{args.chain}'
                                      f'{old_to_new[orig] if args.renumber else orig}'}
                     for aa, orig in checked],
        'model_hotspot_string': hot_model,
        'notes': notes,
        'outputs': {'target_pdb': os.path.abspath(pdb_path),
                    'hotspots_txt': os.path.abspath(hot_path)},
    }

    # --- round trip: read back what was written --------------------------
    rt = read_structure(pdb_path)
    rt_table, rt_ids = residue_table(rt)
    rt_chains = sorted({a['chain'] for a in rt})
    rt_problems = []
    if rt_chains != [args.chain]:
        rt_problems.append(f'written file holds chains {rt_chains}')
    if len(rt_ids) != len(resids):
        rt_problems.append(f'written file holds {len(rt_ids)} residues, expected {len(resids)}')
    if len(rt) != len(kept):
        rt_problems.append(f'written file holds {len(rt)} atoms, expected {len(kept)}')
    for aa, n in emit:
        if n not in rt_table:
            rt_problems.append(f'hotspot {n} missing from the written file')
        elif rt_table[n][1] != aa:
            rt_problems.append(f'hotspot {aa}{n} reads back as {rt_table[n][1]}{n}')
    for a, b in zip(kept[:], rt):
        if abs(a['x'] - b['x']) > 0.001 or abs(a['y'] - b['y']) > 0.001 \
                or abs(a['z'] - b['z']) > 0.001:
            rt_problems.append('coordinates changed on write')
            break
    record['round_trip'] = 'passed' if not rt_problems else rt_problems
    record['outputs']['target_pdb_sha256'] = sha256(pdb_path)

    with open(json_path, 'w') as fh:
        json.dump(record, fh, indent=2)
        fh.write('\n')

    # --- report -----------------------------------------------------------
    print(f'source        {args.structure}  chain {args.chain}  (chains present: '
          f'{",".join(chains)})')
    print(f'window        {lo}-{hi}  ->  kept {len(resids)} residues '
          f'{resids[0]}-{resids[-1]}, {len(kept)} atoms')
    print(f'gaps in it    {gaps if gaps else "none"}')
    if gaps:
        print('              a structure-conditioned model may treat each gap as a '
              'chain break')
    if offset is not None:
        print(f'reference     offset +{offset}  (original + {offset} = reference '
              f'numbering)')
    print(f'numbering     {"renumbered from 1" if args.renumber else "kept as in the source"}')
    print()
    print(f'{"amino acid":>11} {"in file":>8} {"original":>9} {"reference":>10}  model')
    for aa, orig in checked:
        inf = old_to_new[orig] if args.renumber else orig
        ref = f'{aa}{orig + offset}' if offset is not None else '-'
        print(f'{aa + str(orig):>11} {aa + str(inf):>8} {aa + str(orig):>9} {ref:>10}  '
              f'{args.chain}{inf}')
    print()
    print(f'model string  {hot_model}')
    for n in notes:
        print(f'note          {n}')
    print()
    print(f'wrote         {pdb_path}')
    print(f'              {hot_path}')
    print(f'              {json_path}')
    print(f'round trip    {"passed" if not rt_problems else "FAILED"}')
    for p in rt_problems:
        print(f'              ! {p}')
    return 1 if rt_problems else 0


if __name__ == '__main__':
    sys.exit(main())
