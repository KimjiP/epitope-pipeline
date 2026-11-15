#!/usr/bin/env python3
"""Run the whole epitope pipeline on any target, with one command.

Usage:
    python3 run_pipeline.py --uniprot ACC --name NAME [--region LO-HI]
                            [--limit N|all] [--max-resolution A]
                            [--exclude-uniprot ACC,ACC,...] [--domains NAME=LO-HI,...]
                            [--out DIR] [--cache DIR]
    python3 run_pipeline.py --from-lock inputs.lock.json --out DIR

Example, the same six steps run on EGFR:
    python3 run_pipeline.py --uniprot P00533 --name EGFR \
        --exclude-uniprot P04626,P21860,Q15303 --out runs/egfr

What it does, step by step, printed as it goes:

  1. Read what is known: the sequence, the parts and the list of structures, from UniProt.
  2. Collect the published structures that reach the outside of the protein.
  3. Put them in one coordinate frame. Each file's own residue numbering is recovered by
     alignment; relatives of the target are dropped so a receptor partner is not mistaken for
     a binder.
  4. Count distinct binders, not files, and ask where they converge, and whether that is more
     than chance.
  5. Measure the surface.
  6. Write the outputs.

Steps 3 to 5 are atlas.py, convergence.py and surface_props.py. This script fetches the inputs,
picks a reference structure and reports. Every run writes inputs.lock.json and
inputs.snapshot.json; --from-lock replays them.

Stdlib only.
"""
import sys, os, re, json, hashlib, datetime, argparse, subprocess, urllib.request, urllib.error
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

UA = {'User-Agent': 'epitope-pipeline/1.0'}
PDB_ID = re.compile(r'^[0-9][A-Za-z0-9]{3}$')
MAX_AUTO = 100        # more structures than this and the user must say how many to read


def get(url, timeout=60):
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


def banner(n, title, detail=''):
    print()
    print(f'\033[1mSTEP {n}  {title}\033[0m' if sys.stdout.isatty()
          else f'STEP {n}  {title}')
    if detail:
        for line in detail.split('\n'):
            print(f'        {line}')


def uniprot_entry(acc):
    data = json.loads(get(f'https://rest.uniprot.org/uniprotkb/{acc}.json'))
    seq = data['sequence']['value']
    name = (data.get('proteinDescription', {}).get('recommendedName', {})
                .get('fullName', {}).get('value', acc))
    pdbs = []
    for x in data.get('uniProtKBCrossReferences', []):
        if x.get('database') != 'PDB' or not PDB_ID.match(x.get('id', '')):
            continue
        rng = None
        for p in x.get('properties', []):
            if p.get('key') == 'Chains':
                m = re.search(r'=(\d+)-(\d+)', p.get('value', ''))
                if m:
                    rng = (int(m.group(1)), int(m.group(2)))
        pdbs.append({'id': x['id'], 'range': rng})

    # Domain boundaries, in UniProt numbering, from the entry's own annotation.
    # Without these the atlas would label this target's residues with the
    # domain names of whatever target ran last. A feature whose description is
    # only "Approximate" (EGFR has two) carries no name and is left out.
    def bounds(f):
        loc = f.get('location', {})
        a = loc.get('start', {}).get('value')
        b = loc.get('end', {}).get('value')
        return (int(a), int(b)) if a is not None and b is not None else (None, None)

    # EVERY extracellular segment is kept. A membrane protein that crosses the membrane several
    # times (a GPCR, CD20) has several, and the first version kept only the last one.
    doms, seen, sig, regions = [], 0, None, []
    for f in data.get('features', []):
        t, desc = f.get('type'), (f.get('description') or '')
        lo, hi = bounds(f)
        if t == 'Signal' and hi:
            sig = hi
        elif t == 'Topological domain' and 'xtracellular' in desc and lo:
            regions.append((lo, hi))
        elif t in ('Domain', 'Repeat') and lo and desc.strip().lower() != 'approximate':
            seen += 1
            label = re.sub(r'[^A-Za-z0-9]+', '-', desc or f'D{seen}').strip('-')
            doms.append((f'{label[:18] or "D"}-{seen}', lo, hi))
    return seq, name, pdbs, doms, sig, regions


def parse_domains(text):
    """'I=23-217,II=218-341' -> [('I', 23, 217), ('II', 218, 341)]"""
    out = []
    for part in filter(None, (text or '').split(',')):
        name, _, rng = part.partition('=')
        lo, _, hi = rng.partition('-')
        try:
            out.append((name.strip(), int(lo), int(hi)))
        except ValueError:
            sys.exit(f'--domains: cannot read "{part}"; write it as NAME=LO-HI')
    return out


def fmt_regions(regions):
    return ', '.join(f'{lo}-{hi}' for lo, hi in regions)


def in_regions(pdbs, regions):
    """The structures that cover any part of any of the regions (UniProt numbering)."""
    return [p for p in pdbs if p['range'] and any(
        p['range'][1] >= lo and p['range'][0] <= hi for lo, hi in regions)]


def download(pdb_id, cache):
    if not PDB_ID.match(pdb_id):
        return pdb_id, False, 'not a PDB id'
    dest = os.path.join(cache, f'{pdb_id}.cif')
    if os.path.exists(dest) and os.path.getsize(dest) > 2000:
        return pdb_id, True, 'cached'
    try:
        blob = get(f'https://files.rcsb.org/download/{pdb_id}.cif', timeout=120)
    except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as e:
        return pdb_id, False, str(e)[:40]
    if not blob.lstrip().startswith(b'data_'):
        return pdb_id, False, 'not an mmCIF file'
    # write beside the target and rename, so a run that is killed mid-write cannot leave a
    # truncated file that the next run would take for a cached structure
    tmp = f'{dest}.{os.getpid()}.part'
    with open(tmp, 'wb') as fh:
        fh.write(blob)
    os.replace(tmp, dest)
    return pdb_id, True, 'downloaded'


def sha256_of(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


ENTITY_QUERY = '''{
  entries(entry_ids: [%s]) {
    rcsb_id
    struct { title }
    rcsb_entry_info { resolution_combined experimental_method }
    rcsb_accession_info { initial_release_date }
    polymer_entities {
      rcsb_polymer_entity { pdbx_description }
      entity_poly { rcsb_entity_polymer_type }
      rcsb_polymer_entity_container_identifiers { auth_asym_ids entity_id }
      rcsb_entity_source_organism { ncbi_scientific_name }
    }
  }
}'''


def canonical_entities(doc):
    '''Put the PDB's reply in a fixed order.

    The service returns the entities of one entry in a different order from one
    call to the next, so two downloads of the same data differ byte for byte.
    Sorting makes the file, and any hash of it, stable.'''
    entries = (doc.get('data') or {}).get('entries') or []
    for e in entries:
        pes = e.get('polymer_entities') or []
        for pe in pes:
            ids = pe.get('rcsb_polymer_entity_container_identifiers') or {}
            if isinstance(ids.get('auth_asym_ids'), list):
                ids['auth_asym_ids'] = sorted(ids['auth_asym_ids'])
            if isinstance(pe.get('rcsb_entity_source_organism'), list):
                pe['rcsb_entity_source_organism'].sort(
                    key=lambda o: json.dumps(o, sort_keys=True))

        def eid(pe):
            try:
                return int((pe.get('rcsb_polymer_entity_container_identifiers') or {})
                           .get('entity_id'))
            except (TypeError, ValueError):
                return 0
        pes.sort(key=eid)
    entries.sort(key=lambda e: e.get('rcsb_id', ''))
    return doc


def fetch_entities(ids, path):
    """What each chain in each file actually is, from the PDB's own metadata.

    Without this every binder is reported as "unnamed", which makes the output
    unreadable: the point of the atlas is to say *which* binders exist."""
    q = ENTITY_QUERY % ', '.join(f'"{i}"' for i in ids)
    body = json.dumps({'query': q}).encode()
    req = urllib.request.Request('https://data.rcsb.org/graphql', data=body,
                                 headers={**UA, 'Content-Type': 'application/json'})
    with urllib.request.urlopen(req, timeout=120) as r:
        doc = json.loads(r.read())
    if not (doc.get('data') or {}).get('entries'):
        return None, 0
    doc = canonical_entities(doc)
    with open(path, 'w') as fh:
        json.dump(doc, fh, sort_keys=True)
    return path, len(doc['data']['entries'])


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--uniprot', help='accession, e.g. P00533 (read from the lock on a replay)')
    ap.add_argument('--name', help='short name for the target (read from the lock on a replay)')
    ap.add_argument('--region', help='LO-HI in UniProt numbering; keep only '
                                     'structures overlapping it')
    ap.add_argument('--limit', default=None, metavar='N|all',
                    help='how many structures to read: a number, or "all". Default: all, unless '
                         f'more than {MAX_AUTO} reach the outside of the protein, when you must '
                         'choose. A number takes the FIRST N in PDB-id order, which is roughly '
                         'the oldest first, so the counts are then a sample.')
    ap.add_argument('--max-resolution', type=float, default=4.0, metavar='A',
                    help='skip structures solved at a resolution worse than this, in angstrom '
                         '(default 4.0; 0 keeps all). NMR entries have none and are kept.')
    ap.add_argument('--exclude-uniprot', default='',
                    help='relatives of the target, comma separated; their chains '
                         'are receptor partners, not binders')
    ap.add_argument('--domains', metavar='NAME=LO-HI,...',
                    help='name the antigen\'s domains, in UniProt numbering, instead of using '
                         'the ones UniProt annotates (it names none for HER2 or EGFR)')
    ap.add_argument('--from-lock', metavar='LOCK',
                    help='replay a previous run exactly: the same PDB entries and the '
                         'same reference structure, with every file checked against '
                         'the hash recorded then')
    ap.add_argument('--out', default='demo-out')
    ap.add_argument('--cache', help='where structures go (default OUT/pdb-cache)')
    args = ap.parse_args()

    out = os.path.abspath(args.out)
    cache = os.path.abspath(args.cache) if args.cache else os.path.join(out, 'pdb-cache')
    os.makedirs(out, exist_ok=True)      # not just a side effect of making the cache
    os.makedirs(cache, exist_ok=True)

    lock = json.load(open(args.from_lock)) if args.from_lock else None
    snap = None
    if lock:
        if args.uniprot and args.uniprot != lock['uniprot']:
            sys.exit(f'--uniprot {args.uniprot} does not match the lock ({lock["uniprot"]})')
        bad = [x['id'] for x in lock['structures'] if not PDB_ID.match(str(x.get('id')))]
        if bad:
            sys.exit(f'the lock holds ids that are not PDB ids: {bad[:3]}')
        args.uniprot = lock['uniprot']
        args.name = args.name or lock['target']
        print(f'Replaying {args.from_lock}: {len(lock["structures"])} structures '
              f'recorded on {lock["generated"][:10]}')
        sn = lock.get('snapshot')
        if sn:
            sp = os.path.join(os.path.dirname(os.path.abspath(args.from_lock)), sn['file'])
            if not os.path.exists(sp):
                sys.exit(f'the lock names {sn["file"]} next to it, and it is missing')
            if sha256_of(sp) != sn['sha256']:
                sys.exit(f'{sn["file"]} does not match the hash in the lock; it was edited')
            snap = json.load(open(sp))
    elif not (args.uniprot and args.name):
        ap.error('--uniprot and --name are required unless --from-lock is given')

    print(f'Target: {args.name} ({args.uniprot})   output: {out}')

    # ---------------------------------------------------------------- 1 ----
    banner(1, 'Read what is already known about the protein',
           'its sequence, its parts and its list of known structures, '
           'from the UniProt entry')
    if snap:
        u = snap['uniprot']
        seq, fullname, sig = u['sequence'], u['full_name'], u['signal']
        pdbs = [{'id': p['id'], 'range': tuple(p['range']) if p['range'] else None}
                for p in u['pdbs']]
        doms = [(n, a, b) for n, a, b in u['domains']]
        regions = ([tuple(r) for r in u['regions']] if 'regions' in u
                   else ([tuple(u['ecd'])] if u.get('ecd') else []))
        print('        read from the snapshot beside the lock; no network call')
    else:
        seq, fullname, pdbs, doms, sig, regions = uniprot_entry(args.uniprot)
        if args.domains:
            doms = parse_domains(args.domains)
        if lock:
            print('        ! this lock has no snapshot, so UniProt is read live and a change '
                  'there could change the result')
    fasta = os.path.join(out, f'{args.uniprot}.fasta')
    with open(fasta, 'w') as fh:
        fh.write(f'>{args.uniprot}\n')
        for i in range(0, len(seq), 60):
            fh.write(seq[i:i+60] + '\n')
    print(f'        {fullname}')
    print(f'        sequence {len(seq)} residues  ->  {fasta}')
    print(f'        {len(pdbs)} structures listed on the entry')
    seq_sha = hashlib.sha256(seq.encode()).hexdigest()
    if lock and lock['sequence_sha256'] != seq_sha:
        sys.exit('the UniProt sequence for this accession has changed since the lock '
                 'was written; the recorded residue numbers no longer apply')
    if sig:
        print(f'        signal peptide 1-{sig}: files may number from {sig + 1} '
              f'or from 1, so numbering is recovered per file in step 3')
    if doms:
        print(f'        {len(doms)} named domains: '
              + ', '.join(f'{n} {a}-{b}' for n, a, b in doms[:4])
              + (' ...' if len(doms) > 4 else ''))
    else:
        print('        no named domains; the atlas will not name domains '
              '(--domains NAME=LO-HI,... adds them)')

    # ---------------------------------------------------------------- 2 ----
    banner(2, 'Collect the published structures of it',
           'the list on the UniProt entry; download them')
    why = ''
    if lock:
        regions_used = [tuple(r) for r in (lock.get('regions')
                                           or ([lock['region']] if lock.get('region') else []))]
        limit = lock.get('limit')
        max_res = lock.get('max_resolution') or 0.0
    else:
        if args.region:
            regions_used = [tuple(int(v) for v in args.region.split('-'))]
            why = 'asked for'
        else:
            regions_used = regions
            why = 'the outside-the-cell part, per UniProt'
        limit = None if args.limit is None or str(args.limit).lower() == 'all' else int(args.limit)
        max_res = args.max_resolution or 0.0

    chosen = pdbs
    if lock:
        byid = {p['id']: p for p in pdbs}
        chosen = [byid.get(x['id'], {'id': x['id'], 'range': None}) for x in lock['structures']]
        print(f'        the lock names {len(chosen)} structures; the entry listed '
              f'{len(pdbs)}')
        gone = [x['id'] for x in lock['structures'] if x['id'] not in byid]
        if gone:
            print('        ! no longer on the entry, fetching anyway: ' + ', '.join(gone))
        in_region = lock.get('structures_in_region') or len(chosen)
    else:
        if regions_used:
            chosen = in_regions(pdbs, regions_used)
            print(f'        keeping the {len(chosen)} of {len(pdbs)} that overlap '
                  f'residues {fmt_regions(regions_used)} ({why})')
            print('        an antibody drug must reach its target from outside the cell')
        in_region = len(chosen)
        if args.limit is None and in_region > MAX_AUTO:
            sys.exit(f'{in_region} structures reach the outside of the protein, and reading all '
                     f'of them is a large download. Run again with --limit all to read them '
                     f'all, or --limit N to read the first N (oldest first).')
        if limit is not None and in_region > limit:
            print(f'        taking the first {limit} of {in_region}, in PDB-id order (roughly '
                  f'the oldest first). The rest are NOT read: the counts below are a sample.')
            chosen = chosen[:limit]
    ids = [p['id'] for p in chosen]
    ok, failed = [], []
    with ThreadPoolExecutor(max_workers=8) as ex:
        for pid, good, how in ex.map(lambda p: download(p, cache), ids):
            if good:
                ok.append(pid)
            else:
                failed.append(pid)
                print(f'        ! {pid}: {how}')
    print(f'        {len(ok)} structure files in {cache}')
    if failed:
        print(f'        ! {len(failed)} structures could not be fetched and are left out: '
              + ', '.join(failed))
    if len(ok) < 3:
        sys.exit('not enough structures to continue')
    # The hash of every file this run actually used. Written to the lock, and
    # checked against it on a replay: a structure can be re-deposited, and then
    # the same PDB id is no longer the same coordinates.
    hashes = {pid: sha256_of(os.path.join(cache, f'{pid}.cif')) for pid in ok}
    if lock:
        was = {x['id']: x['sha256'] for x in lock['structures']}
        moved = [p for p in ok if p in was and was[p] != hashes[p]]
        if moved:
            print('        ! re-deposited since the lock, contents differ: '
                  + ', '.join(moved))
        else:
            print(f'        every file matches the hash in the lock')

    excl_path = None
    # On a replay the lock is authoritative. Leaving the exclusions out brought
    # six receptor-partner entries back in as binders and moved the answer.
    accs = ([a.strip() for a in args.exclude_uniprot.split(',') if a.strip()]
            if not lock else list(lock.get('exclude_uniprot') or []))
    if lock and accs:
        print(f'        relatives from the lock: {", ".join(accs)}')
    rel_seqs = {}
    if accs:
        excl_path = os.path.join(out, 'relatives.fasta')
        for a in accs:
            rel_seqs[a] = (snap['relatives'][a] if snap and a in snap.get('relatives', {})
                           else get(f'https://rest.uniprot.org/uniprotkb/{a}.fasta').decode())
        with open(excl_path, 'w') as fh:
            for a in accs:
                fh.write(rel_seqs[a])
        print(f'        relatives to exclude: {", ".join(accs)}')

    # ---------------------------------------------------------------- 3 ----
    banner(3, 'Put them all in one coordinate frame',
           "recover each file's own residue numbering, then superpose every "
           'complex')
    from atlas import parse, chain_residues, intkeyed, find_offset, AA3   # noqa: E402
    import convergence                                                     # noqa: E402
    best = None
    for pid in ok:
        try:
            ch = chain_residues(parse(os.path.join(cache, f'{pid}.cif')))
        except Exception:
            continue
        for c, res in ch.items():
            ir = intkeyed(res)
            if len(ir) < 50:
                continue
            o = find_offset(ir, seq)
            if o is None:
                continue
            if best is None or len(ir) > best[2]:
                best = (pid, c, len(ir), o)
    if not best:
        sys.exit('no chain in any file matches the reference sequence')
    ref_pdb, ref_chain, ref_n, ref_off = best
    if lock:
        lr = lock['reference']
        if (lr['pdb'], lr['chain']) != (ref_pdb, ref_chain):
            print(f'        using the locked reference {lr["pdb"]} chain {lr["chain"]} '
                  f'(today the scan would pick {ref_pdb} chain {ref_chain})')
            ref_pdb, ref_chain, ref_off = lr['pdb'], lr['chain'], lr['offset']
            ref_n = lr.get('modelled', ref_n)
    print(f'        reference: {ref_pdb} chain {ref_chain}, {ref_n} modelled '
          f'residues, offset {ref_off:+d}')

    # atlas.py takes domain bounds in the REFERENCE chain's own numbering, and
    # its default is HER2's. uniprot_index = resi + ref_off, so resi = i - off.
    dom_arg = ','.join(f'{n}={a - ref_off}-{b - ref_off}'
                       for n, a, b in doms if b - ref_off > 0)
    atlas_out = os.path.join(out, 'atlas-out')
    ref_cif = os.path.join(cache, f'{ref_pdb}.cif')
    cmd = [sys.executable, '-I', os.path.join(HERE, 'atlas.py'),
           ref_cif, ref_chain, fasta, cache,
           '--domains', dom_arg, '--only', ','.join(ok), '--out', atlas_out,
           '--max-resolution', f'{max_res:g}']
    if excl_path:
        cmd += ['--exclude', excl_path]
    ent_doc = None
    try:
        if snap and snap.get('entities'):
            ent_path = os.path.join(out, 'entities.json')
            with open(ent_path, 'w') as fh:
                json.dump(snap['entities'], fh, sort_keys=True)
            n_ent = len(snap['entities']['data']['entries'])
        else:
            ent_path, n_ent = fetch_entities(ok, os.path.join(out, 'entities.json'))
    except Exception as e:
        ent_path, n_ent = None, 0
        print(f'        ! chain names unavailable ({str(e)[:50]}); '
              f'binders will be unnamed')
    if ent_path:
        ent_doc = json.load(open(ent_path))
        cmd += ['--entities', ent_path]
        print(f'        chain names and resolutions for {n_ent} entries from '
              f'{"the snapshot" if snap else "the PDB"}')
    elif max_res:
        print('        ! no resolutions available, so --max-resolution cannot be applied')
    print(f'        running atlas.py over {len(ok)} files ...')
    r = subprocess.run(cmd, capture_output=True, text=True)
    atlas_log = os.path.join(out, 'atlas.log')
    with open(atlas_log, 'w') as fh:
        fh.write(r.stdout + '\n' + r.stderr)
    if r.returncode != 0:
        print(r.stdout[-2500:]); print(r.stderr[-2500:])
        sys.exit('atlas.py failed')
    for line in r.stdout.splitlines():
        if 'entries used' in line or 'excluding chains' in line:
            print(f'        {line.strip()}')

    A = json.load(open(os.path.join(atlas_out, 'atlas.json')))
    entries = [e for e in A.get('entries', []) if e.get('status') == 'ok']
    conv = {}
    for e in entries:
        conv[e.get('numbering', '?')] = conv.get(e.get('numbering', '?'), 0) + 1
    print(f'        {len(entries)} complexes placed in one frame')
    for k, v in sorted(conv.items(), key=lambda kv: -kv[1]):
        print(f'          {v:>3}  {k}')
    n_lowres = sum(1 for e in A.get('entries', [])
                   if str(e.get('status', '')).startswith('resolution worse'))
    if n_lowres:
        print(f'        {n_lowres} entries skipped: solved at a resolution worse than '
              f'{max_res:g} A (--max-resolution 0 keeps them)')

    # ---------------------------------------------------------------- 4 ----
    banner(4, 'Count distinct binders, not files',
           'one antibody appears in six entries; collapse it to one')
    nd = A.get('n_distinct')
    nu = A.get('n_units')
    units = A.get('units', [])
    clusters = A.get('clusters', [])
    kinds = {k: sum(1 for c in clusters if c.get('kind') == k)
             for k in ('antibody', 'scaffold', 'other')}
    print(f'        {len(entries)} complexes  ->  {nu} binding units  ->  '
          f'{nd} distinct binders')
    print(f'        of the {nd}: {kinds["antibody"]} antibodies and nanobodies, '
          f'{kinds["scaffold"]} engineered scaffolds, {kinds["other"]} other proteins')
    print('        (one binder solved six times is one binder, not six)')
    fmts = {}
    for u in units:
        f = u.get('format') or '?'
        fmts[f] = fmts.get(f, 0) + 1
    print('        shapes of the binding units:')
    for k, v in sorted(fmts.items(), key=lambda kv: -kv[1])[:6]:
        print(f'          {v:>3}  {k}')
    # Not everything that binds is a drug. Other proteins are often the target's own ligands or
    # receptors. They are counted and listed, and kept out of the patch below, which asks where
    # ANTIBODIES converge.
    other = [c for c in clusters if c.get('kind') == 'other']
    if other:
        print(f'        {len(other)} distinct binders are other proteins (often the antigen\'s '
              f'own ligands or receptors):')
        for c in other[:6]:
            print(f'          {(c.get("display") or c["label"])[:46]}')
    repeats = sorted((c for c in clusters if len(c.get('members', [])) > 1),
                     key=lambda c: -len(c['members']))
    if repeats:
        print('        solved more than once:')
        for c in repeats[:5]:
            print(f'          {len(c["members"]):>2} x  '
                  f'{(c.get("display") or c.get("label") or c["id"])[:34]:<34} '
                  f'{", ".join(c.get("pdbs", [])[:5])}')
    n_misfit = sum(1 for u in units if u.get('rmsd', 0) > 3.0)
    n_quat = sum(1 for u in units if u.get('other_copy_share', 0) >= 0.2)
    if n_misfit:
        print(f'        ! {n_misfit} units are placed with a fit error above 3 A; their position '
              f'in the 3D page is approximate')
    if n_quat:
        print(f'        ! {n_quat} units make 20% or more of their contacts with ANOTHER copy of '
              f'the antigen. Only one copy is read, so their epitope is cut short.')
    # Antigen residues with an insertion code (52A beside 52) are not read; HIV Env and some
    # antibody-numbered antigens use them. None in HER2 or EGFR.
    n_ic = sum(e.get('antigen_residues_with_insertion_codes_skipped', 0) for e in entries)
    if n_ic:
        print(f'        ! {n_ic} antigen residues carry insertion codes and are skipped, so any '
              f'contact on them is missed')

    # The surface is measured first because the baseline below needs to know which residues are
    # exposed. Its lines are printed under step 5.
    props = os.path.join(out, 'props.json')
    r2 = subprocess.run([sys.executable, '-I', os.path.join(HERE, 'surface_props.py'),
                         ref_cif, ref_chain,
                         '--atlas', os.path.join(atlas_out, 'atlas.json'),
                         '--json', props],
                        capture_output=True, text=True)
    rsasa = {}
    if r2.returncode == 0 and os.path.exists(props):
        rsasa = {int(k): v['rsasa'] for k, v in json.load(open(props))['residues'].items()}
    exposed = {r for r, x in rsasa.items() if x >= 0.10}

    # Where do the binders converge? Footprints are the residues a binder touches in at least
    # half of its structures. See convergence.py for the rule and the baseline.
    refres = intkeyed(chain_residues(parse(ref_cif))[ref_chain])
    ca = {r: ats['CA'] for r, (comp, ats) in refres.items() if 'CA' in ats}
    dens = {int(k): v for k, v in A.get('density_distinct', {}).items()}
    foot = {c['id']: set(c.get('epitope_core') or []) for c in clusters}
    names = {c['id']: (c.get('display') or c.get('label') or c['id']) for c in clusters}
    nb = convergence.neighbours(ca)

    def dom_name(uniprot_number):
        for n, a, b in doms:
            if a <= uniprot_number <= b:
                return n
        return None

    patches = {}
    for basis, bids in (('antibodies and nanobodies', {c['id'] for c in clusters
                                                       if c.get('kind') == 'antibody'}),
                        ('all binders', {c['id'] for c in clusters})):
        fp = {i: foot[i] for i in bids}
        res = convergence.find_patch(fp, ca, nb)
        if res and exposed:
            res['chance'] = convergence.chance(
                convergence.baseline(fp, ca, exposed, nb), res['hits'])
        patches[basis] = res
    primary = 'antibodies and nanobodies' if kinds['antibody'] >= 4 else 'all binders'
    hot = patches[primary]

    def say_patch(basis, res):
        if not res:
            print(f'        {basis}: too few binders with 3 or more contact residues to look '
                  f'for a patch')
            return
        c = res['centre'] + ref_off
        d = dom_name(c)
        print(f'        most-targeted patch, {basis}: {res["hits"]} of {res["considered"]} '
              f'touch the sphere (radius {convergence.RADIUS:g} A) around residue {c} '
              f'(UniProt)' + (f', {d}' if d else ''))
        ch_ = res.get('chance')
        if ch_:
            verdict = ('above chance' if ch_['p'] < 0.05 else
                       'NOT above chance: these binders are spread over the antigen')
            print(f'          chance: random footprints of the same sizes give a median of '
                  f'{ch_["median"]} (95th percentile {ch_["p95"]}); p = {ch_["p"]:g}; {verdict}')
        if basis == primary:
            print(f'          {", ".join(names[i][:22] for i in res["binders"][:6])}'
                  + (' ...' if len(res['binders']) > 6 else ''))
    say_patch(primary, hot)
    if primary != 'all binders':
        say_patch('all binders', patches['all binders'])

    # ---------------------------------------------------------------- 5 ----
    banner(5, 'Measure the surface',
           'exposure, hydrophobicity, burial, and how many binders touch each amino acid')
    if r2.returncode != 0:
        print(r2.stdout[-1500:]); print(r2.stderr[-1500:])
    else:
        for line in r2.stdout.splitlines():
            if any(k in line for k in ('exposed', 'hydrophobic', 'burial', 'rsasa',
                                       'Bound (', 'property')):
                print(f'        {line.strip()}')

    # ---------------------------------------------------------------- 6 ----
    banner(6, 'Hand over what a design model needs',
           'a hotspot file a design model can read, and the measured data behind it')
    hot_files, spec, page = {}, '', None
    supported = bool(hot and hot.get('chance') and hot['chance']['p'] < 0.05)
    if hot:
        cdb = defaultdict(int)
        for i in [c['id'] for c in clusters if (c.get('kind') == 'antibody')
                  or primary == 'all binders']:
            for r_ in foot[i]:
                cdb[r_] += 1
        picked = convergence.pick_hotspots(hot['touched'], cdb, rsasa, dens)
        if not picked:
            print('        no residue in the patch is touched by 2 or more binders and at '
                  'least 20% exposed, so no hotspot file is written')
        else:
            spec = ','.join(f'{AA3.get(refres[r_][0], "X")}{r_}' for r_ in picked if r_ in refres)
            first, last = min(refres), max(refres)
            lo_w, hi_w = max(first, min(picked) - 20), min(last, max(picked) + 20)
            tdir = os.path.join(out, 'target')
            stem = args.name.lower()
            r3 = subprocess.run(
                [sys.executable, '-I', os.path.join(HERE, 'prepare_target.py'),
                 ref_cif, ref_chain,
                 '--window', f'{lo_w}-{hi_w}',
                 '--hotspots', spec, '--uniprot', fasta,
                 '--outdir', tdir, '--name', stem],
                capture_output=True, text=True)
            if r3.returncode == 0:
                print(f'        hotspots, in {ref_pdb} chain {ref_chain} numbering: {spec}')
                print(f'          the same residues in UniProt numbering: '
                      + ','.join(f'{AA3.get(refres[r_][0], "X")}{r_ + ref_off}'
                                 for r_ in picked if r_ in refres))
                print(f'          (this file numbers from {ref_off:+d} against '
                      f'UniProt; the design model needs the file\'s own numbers)')
                print('          each is touched by 2 or more binders and is at least 20% '
                      'exposed in the free antigen')
                if not supported:
                    print('          ! the patch is NOT above chance (step 4), so these residues are a '
                          'starting guess, not a consensus')
                for f in sorted(os.listdir(tdir)):
                    if f.startswith(stem) and '_full' not in f:
                        hot_files[f] = os.path.join(tdir, f)
                        print(f'          {os.path.join(tdir, f)}')
            else:
                print('        ! prepare_target.py refused the hotspot list, '
                      'which is what it is for:')
                for line in (r3.stdout + r3.stderr).strip().splitlines()[-4:]:
                    print(f'          {line.strip()}')

        # The same 3D page, with this target's numbers in it.
        tdir = os.path.join(out, 'target')
        stem = args.name.lower()
        r4 = subprocess.run(
            [sys.executable, '-I', os.path.join(HERE, 'prepare_target.py'),
             ref_cif, ref_chain,
             '--window', f'{min(refres)}-{max(refres)}',
             '--hotspots', spec or f'{AA3.get(refres[min(refres)][0], "X")}{min(refres)}',
             '--uniprot', fasta, '--outdir', tdir, '--name', f'{stem}_full'],
            capture_output=True, text=True)
        whole = os.path.join(tdir, f'{stem}_full_target.pdb')
        if r4.returncode == 0 and os.path.exists(whole):
            import demo_page                                          # noqa: E402
            cdist = {k: v for k, v in A.get('counts_distinct', {}).items()}
            c0 = hot['centre'] + ref_off
            d0 = dom_name(c0)
            page = demo_page.write(
                os.path.join(out, 'viewer'),
                {'title': f'{args.name}: where binders land',
                 'subtitle': f'{len(ok)} published structures, {nd} distinct '
                             f'binders, one coordinate frame. '
                             f'Reference {ref_pdb} chain {ref_chain}.',
                 'density': cdist,
                 'facts': [['structures listed', len(pdbs)],
                           ['structures used', len(entries)],
                           ['binding units', nu],
                           ['distinct binders', nd],
                           ['of them antibodies', kinds['antibody']],
                           ['residues measured', len(refres)],
                           ['numbering vs UniProt', f'{ref_off:+d}']],
                 'patch': {'label': f'{hot["hits"]} of {hot["considered"]} {primary} touch '
                                    f'the region around residue {c0} (UniProt)'
                                    + (f', {d0}' if d0 else ''),
                           'resi': hot['residues'],
                           'binders': [names[i] for i in hot['binders'][:8]]}},
                whole, os.path.join(atlas_out, 'binders.pdb'))
            print(f'        3D page: {page}')

    print('        the measured data:')

    def patch_record(basis, res):
        if not res:
            return None
        c = res['centre'] + ref_off
        touched = [r_ + ref_off for r_ in res['touched']]
        return {'basis': basis, 'centre_uniprot': c, 'radius_A': convergence.RADIUS,
                'domain': dom_name(c), 'hits': res['hits'], 'of': res['considered'],
                'rule': f'a binder hits when at least {convergence.MIN_RES} of its contact '
                        f'residues, and at least {int(convergence.MIN_SHARE * 100)}% of them, '
                        f'lie in the sphere',
                'binders': [names[i] for i in res['binders']],
                'touched_uniprot': touched,
                'touched_runs_uniprot': convergence.runs_text(touched),
                'chance': res.get('chance')}
    released = sorted(e.get('rcsb_accession_info', {}).get('initial_release_date', '')[:10]
                      for e in ((ent_doc or {}).get('data') or {}).get('entries') or []
                      if e['rcsb_id'] in set(ok) and e.get('rcsb_accession_info'))
    summary = {'target': args.name, 'uniprot': args.uniprot, 'full_name': fullname,
               'sequence_length': len(seq), 'signal_peptide': sig,
               'structures_listed': len(pdbs), 'structures_in_region': in_region,
               'structures_used': len(ok),
               'sampled': bool(limit is not None and in_region > len(ok) + len(failed)),
               'limit': limit, 'max_resolution': max_res or None,
               'released_through': released[-1] if released else None,
               'reference': f'{ref_pdb} chain {ref_chain}',
               'reference_offset': ref_off,
               'complexes': len(entries), 'units': nu, 'distinct_binders': nd,
               'distinct_antibodies': kinds['antibody'],
               'distinct_scaffolds': kinds['scaffold'], 'distinct_other': kinds['other'],
               'entries_skipped': A.get('skipped'),
               'numbering_conventions': conv, 'unit_formats': fmts,
               'domains_uniprot': [{'name': n, 'first': a, 'last': b}
                                   for n, a, b in doms],
               'most_targeted_patch': patch_record(primary, hot),
               'patch_by_basis': {b: patch_record(b, p) for b, p in patches.items()},
               'quality': {'units_fit_error_over_3A': n_misfit,
                           'units_with_contacts_on_another_copy': n_quat,
                           'antigen_residues_skipped_for_insertion_codes': n_ic,
                           'downloads_failed': failed},
               'hotspots_uniprot': ([f'{AA3.get(refres[r_][0], "X")}{r_ + ref_off}'
                                     for r_ in sorted(int(x[1:]) for x in spec.split(','))
                                     if r_ in refres] if spec else []),
               'hotspot_rule': 'touched by 2 or more binders and at least 20% exposed',
               'hotspots_supported': supported,
               'hotspot_files': hot_files, 'viewer': page,
               'outputs': {'atlas': os.path.join(atlas_out, 'atlas.json'),
                           'binders': os.path.join(atlas_out, 'binders.pdb'),
                           'props': props, 'sequence': fasta,
                           'atlas_log': atlas_log}}
    if os.path.exists(props):
        p = json.load(open(props))
        summary['correlations'] = p.get('correlations')
    # The lock is what makes a rerun the same run: the exact entries, the hash of
    # each file, the reference that was chosen, and the sequence they were
    # numbered against. Replay it with --from-lock.
    # The snapshot holds everything else a run reads from the network: the UniProt
    # facts, the sequences of the relatives, and the PDB's chain names. With it a
    # replay needs the structure files and nothing more.
    snap_out = os.path.join(out, 'inputs.snapshot.json')
    with open(snap_out, 'w') as fh:
        json.dump({'schema': 'epitope-pipeline-snapshot/2',
                   'uniprot': {'accession': args.uniprot, 'full_name': fullname,
                               'sequence': seq, 'signal': sig,
                               'domains': [list(d) for d in doms],
                               'regions': [list(r_) for r_ in regions],
                               'pdbs': [{'id': p['id'],
                                         'range': list(p['range']) if p['range'] else None}
                                        for p in pdbs]},
                   'relatives': rel_seqs,
                   'entities': ent_doc},
                  fh, sort_keys=True, separators=(',', ':'))
        fh.write('\n')

    lock_out = os.path.join(out, 'inputs.lock.json')
    with open(lock_out, 'w') as fh:
        json.dump({'schema': 'epitope-pipeline-lock/3',
                   'generated': datetime.datetime.now(datetime.timezone.utc)
                                .replace(microsecond=0).isoformat(),
                   'target': args.name, 'uniprot': args.uniprot,
                   'sequence_sha256': seq_sha, 'sequence_length': len(seq),
                   'regions': [list(r_) for r_ in regions_used],
                   'limit': limit, 'structures_in_region': in_region,
                   'max_resolution': max_res or None,
                   'exclude_uniprot': accs,
                   'reference': {'pdb': ref_pdb, 'chain': ref_chain,
                                 'offset': ref_off, 'modelled': ref_n},
                   'structures': [{'id': p, 'sha256': hashes[p]} for p in ok],
                   'snapshot': {'file': os.path.basename(snap_out),
                                'sha256': sha256_of(snap_out)}},
                  fh, indent=2)
        fh.write('\n')
    summary['lock'] = lock_out
    summary['snapshot'] = snap_out

    # Paths relative to the output directory. An absolute path records whose
    # laptop the run happened on, and these files get committed and published.
    def rel(x):
        return os.path.relpath(x, out) if isinstance(x, str) and os.path.isabs(x) else x
    summary['lock'] = rel(summary['lock'])
    summary['snapshot'] = rel(summary['snapshot'])
    summary['outputs'] = {k: rel(v) for k, v in summary['outputs'].items()}
    summary['hotspot_files'] = {k: rel(v) for k, v in summary['hotspot_files'].items()}
    summary['viewer'] = rel(summary['viewer'])
    print(f'        lock      {lock_out}')

    with open(os.path.join(out, 'summary.json'), 'w') as fh:
        json.dump(summary, fh, indent=2)
        fh.write('\n')
    for k, v in summary['outputs'].items():
        print(f'        {k:<9} {v}')
    print(f'        summary   {os.path.join(out, "summary.json")}')
    print()
    print(f'Same six steps, different target. {args.name}: {len(ok)} structures in, '
          f'{nd} distinct binders out ({kinds["antibody"]} antibodies and nanobodies).')
    return 0


if __name__ == '__main__':
    sys.exit(main())
