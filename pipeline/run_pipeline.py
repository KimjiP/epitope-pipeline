#!/usr/bin/env python3
"""Run the whole epitope pipeline on any target, with one command.

Usage:
    python3 run_pipeline.py --uniprot ACC --name NAME [--region LO-HI]
                            [--limit N] [--exclude-uniprot ACC,ACC,...]
                            [--out DIR] [--cache DIR]

Example — the same six steps the HER2 map shows, run on EGFR:
    python3 run_pipeline.py --uniprot P00533 --name EGFR \
        --region 1-640 --limit 30 \
        --exclude-uniprot P04626,P21860,Q15303 --out ../demo-egfr

Why this exists
---------------
The HER2 analysis was not a one-off. Everything in it is driven by two
identifiers: a UniProt accession and a list of structures. Give it a different
accession and the same steps run. This script is the proof.

What it does, step by step, printed as it goes:

  1. Ask what is already known.   In the real HER2 run this was nine Amass
     queries, and the list of structures came from one field of an Amass gene
     record. Here the same list is read from the UniProt entry's own
     cross-references, so the demo needs no account and no credits.
  2. Collect every published structure.
  3. Put them in one coordinate frame.  Each file's own residue numbering is
     recovered by alignment; relatives of the target are dropped so a receptor
     partner is not mistaken for a binder.
  4. Count distinct binders, not files.
  5. Measure the surface.
  6. Write the outputs.

Steps 3 to 5 are atlas.py and surface_props.py, unchanged. This script only
fetches the inputs, picks a reference structure and reports.

Stdlib only.
"""
import sys, os, re, json, hashlib, datetime, argparse, subprocess, urllib.request, urllib.error
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

UA = {'User-Agent': 'epitope-pipeline/1.0'}


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
        if x.get('database') != 'PDB':
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
    # domain names of whatever target ran last.
    def bounds(f):
        loc = f.get('location', {})
        a = loc.get('start', {}).get('value')
        b = loc.get('end', {}).get('value')
        return (int(a), int(b)) if a is not None and b is not None else (None, None)

    doms, seen, sig, ecd = [], 0, None, None
    for f in data.get('features', []):
        t, desc = f.get('type'), (f.get('description') or '')
        lo, hi = bounds(f)
        if t == 'Signal' and hi:
            sig = hi
        elif t == 'Topological domain' and 'xtracellular' in desc and lo:
            ecd = (lo, hi)
        elif t in ('Domain', 'Repeat') and lo:
            seen += 1
            label = re.sub(r'[^A-Za-z0-9]+', '-', desc or f'D{seen}').strip('-')
            doms.append((f'{label[:18] or "D"}-{seen}', lo, hi))
    return seq, name, pdbs, doms, sig, ecd


def download(pdb_id, cache):
    dest = os.path.join(cache, f'{pdb_id}.cif')
    if os.path.exists(dest) and os.path.getsize(dest) > 2000:
        return pdb_id, True, 'cached'
    try:
        blob = get(f'https://files.rcsb.org/download/{pdb_id}.cif', timeout=120)
    except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as e:
        return pdb_id, False, str(e)[:40]
    with open(dest, 'wb') as fh:
        fh.write(blob)
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
    ap.add_argument('--limit', type=int, default=30,
                    help='how many structures to download (default 30)')
    ap.add_argument('--exclude-uniprot', default='',
                    help='relatives of the target, comma separated; their chains '
                         'are receptor partners, not binders')
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
    banner(1, 'Ask the literature what is already known',
           'For HER2 this was nine Amass queries, and the structure list came\n'
           'from one field of an Amass gene record. Here the same list is read\n'
           'from the UniProt entry, so the demo needs no account.')
    if snap:
        u = snap['uniprot']
        seq, fullname, sig = u['sequence'], u['full_name'], u['signal']
        pdbs = [{'id': p['id'], 'range': tuple(p['range']) if p['range'] else None}
                for p in u['pdbs']]
        doms = [(n, a, b) for n, a, b in u['domains']]
        ecd = tuple(u['ecd']) if u['ecd'] else None
        print('        read from the snapshot beside the lock; no network call')
    else:
        seq, fullname, pdbs, doms, sig, ecd = uniprot_entry(args.uniprot)
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
        print(f'        {len(doms)} annotated domains: '
              + ', '.join(f'{n} {a}-{b}' for n, a, b in doms[:4])
              + (' ...' if len(doms) > 4 else ''))
    else:
        print('        no annotated domains on the entry; the atlas will not '
              'name domains')

    # ---------------------------------------------------------------- 2 ----
    banner(2, 'Collect every published structure of it',
           'one field of one record listed them all; download the lot')
    chosen, region = pdbs, None
    if lock:
        byid = {p['id']: p for p in pdbs}
        chosen = [byid.get(x['id'], {'id': x['id'], 'range': None})
                  for x in lock['structures']]
        region = tuple(lock['region']) if lock.get('region') else None
        print(f'        the lock names {len(chosen)} structures; the entry listed '
              f'{len(pdbs)}')
        gone = [x['id'] for x in lock['structures'] if x['id'] not in byid]
        if gone:
            print('        ! no longer on the entry, fetching anyway: ' + ', '.join(gone))
    elif args.region:
        region = tuple(int(v) for v in args.region.split('-'))
        why = 'asked for'
    elif ecd:
        region = ecd
        why = 'the extracellular part, per UniProt'
    if region and not lock:
        lo, hi = region
        chosen = [p for p in pdbs
                  if p['range'] and p['range'][1] >= lo and p['range'][0] <= hi]
        print(f'        keeping the {len(chosen)} of {len(pdbs)} that overlap '
              f'residues {lo}-{hi} ({why})')
        print('        an antibody drug must reach its target from outside the cell')
    if len(chosen) > args.limit and not lock:
        print(f'        taking the first {args.limit} (raise --limit for more)')
        chosen = chosen[:args.limit]
    ids = [p['id'] for p in chosen]
    ok = []
    with ThreadPoolExecutor(max_workers=8) as ex:
        for pid, good, how in ex.map(lambda p: download(p, cache), ids):
            if good:
                ok.append(pid)
            else:
                print(f'        ! {pid}: {how}')
    print(f'        {len(ok)} structure files in {cache}')
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
    from atlas import parse, chain_residues, intkeyed, find_offset   # noqa: E402
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
    cmd = [sys.executable, '-I', os.path.join(HERE, 'atlas.py'),
           os.path.join(cache, f'{ref_pdb}.cif'), ref_chain, fasta, cache,
           '--domains', dom_arg, '--only', ','.join(ok), '--out', atlas_out]
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

    # ---------------------------------------------------------------- 4 ----
    banner(4, 'Count distinct binders, not files',
           'one antibody appears in six entries; collapse it to one')
    nd = A.get('n_distinct')
    nu = A.get('n_units')
    units = A.get('units', [])
    clusters = A.get('clusters', [])
    print(f'        {len(entries)} complexes  ->  {nu} binding units  ->  '
          f'{nd} distinct binders')
    print('        (one binder solved six times is one binder, not six)')
    fmts = {}
    for u in units:
        f = u.get('format') or '?'
        fmts[f] = fmts.get(f, 0) + 1
    print('        shapes of the binding units:')
    for k, v in sorted(fmts.items(), key=lambda kv: -kv[1])[:6]:
        print(f'          {v:>3}  {k}')
    # Not everything that binds is a drug. On HER2 the receptor partners had to
    # be excluded by sequence; here the target's own growth-factor ligands come
    # through as binders. The pipeline reports them, it does not hide them.
    # A nanobody IS an antibody, so atlas.ANTIBODY_LIKE (which is about
    # paired-chain binders) is too narrow for this question.
    AB = ('Fab', 'scFv', 'antibody fragment', 'VHH', 'single domain (VHH-sized)')
    other = [u for u in units if (u.get('format') or '') not in AB]
    if other:
        print(f'        {len(other)} of {nu} units are not antibody-shaped; '
              f'check whether they belong:')
        for u in other[:6]:
            print(f'          {(u.get("name") or "unnamed")[:40]:<40} '
                  f'{u.get("format")}  ({u["pdb"]})')
    repeats = sorted((c for c in clusters if len(c.get('members', [])) > 1),
                     key=lambda c: -len(c['members']))
    if repeats:
        print('        solved more than once:')
        for c in repeats[:5]:
            print(f'          {len(c["members"]):>2} x  '
                  f'{(c.get("display") or c.get("label") or c["id"])[:34]:<34} '
                  f'{", ".join(c.get("pdbs", [])[:5])}')

    # The most-targeted patch: slide a window and count DISTINCT binders in it,
    # not contacts, so a popular antibody cannot vote twice.
    ep_by_unit = {u['id']: set(u.get('epitope') or []) for u in units}
    cl_ep = {}
    for c in clusters:
        s = set()
        for m in c.get('members', []):
            s |= ep_by_unit.get(m, set())
        if s:
            cl_ep[c.get('display') or c.get('label') or c['id']] = s
    hot = None
    if cl_ep:
        allres = sorted(set().union(*cl_ep.values()))
        W = 24
        for s0 in range(allres[0], allres[-1] - W + 2):
            hits = [n for n, s in cl_ep.items()
                    if sum(1 for r in s if s0 <= r < s0 + W) >= 2]
            if hot is None or len(hits) > hot[2]:
                hot = (s0, s0 + W - 1, len(hits), hits)
    if hot:
        lo, hi, n, who = hot
        print(f'        most-targeted patch: residues {lo + ref_off}-{hi + ref_off} '
              f'(UniProt), hit by {n} of {nd} distinct binders')
        print(f'          {", ".join(w[:22] for w in who[:6])}'
              + (' ...' if len(who) > 6 else ''))

    # ---------------------------------------------------------------- 5 ----
    banner(5, 'Measure the surface and score every site',
           'contacts, exposure, hydrophobicity, and what has precedent')
    props = os.path.join(out, 'props.json')
    r2 = subprocess.run([sys.executable, '-I', os.path.join(HERE, 'surface_props.py'),
                         os.path.join(cache, f'{ref_pdb}.cif'), ref_chain,
                         '--atlas', os.path.join(atlas_out, 'atlas.json'),
                         '--json', props],
                        capture_output=True, text=True)
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
    hot_files = {}
    if hot:
        dens = {int(k): v for k, v in A.get('density_distinct', {}).items()}
        lo, hi = hot[0], hot[1]
        ranked = sorted((r for r in range(lo, hi + 1) if dens.get(r)),
                        key=lambda r: -dens[r])[:8]
        refres = intkeyed(chain_residues(
            parse(os.path.join(cache, f'{ref_pdb}.cif')))[ref_chain])
        from atlas import AA3                                     # noqa: E402
        spec = ','.join(f'{AA3.get(refres[r][0], "X")}{r}'
                        for r in sorted(ranked) if r in refres)
        tdir = os.path.join(out, 'target')
        stem = args.name.lower()
        r3 = subprocess.run(
            [sys.executable, '-I', os.path.join(HERE, 'prepare_target.py'),
             os.path.join(cache, f'{ref_pdb}.cif'), ref_chain,
             '--window', f'{max(1, lo - 20)}-{hi + 20}',
             '--hotspots', spec, '--uniprot', fasta,
             '--outdir', tdir, '--name', stem],
            capture_output=True, text=True)
        if r3.returncode == 0:
            print(f'        hotspots, in {ref_pdb} chain {ref_chain} numbering: '
                  f'{spec}')
            print(f'          the same residues in UniProt numbering: '
                  + ','.join(f'{AA3.get(refres[r][0], "X")}{r + ref_off}'
                             for r in sorted(ranked) if r in refres))
            print(f'          (this file numbers from {ref_off:+d} against '
                  f'UniProt; the design model needs the file\'s own numbers)')
            for f in sorted(os.listdir(tdir)):
                if f.startswith(stem) and '_full' not in f:
                    hot_files[f] = os.path.join(tdir, f)
                    print(f'          {os.path.join(tdir, f)}')
        else:
            print('        ! prepare_target.py refused the hotspot list, '
                  'which is what it is for:')
            for line in (r3.stdout + r3.stderr).strip().splitlines()[-4:]:
                print(f'          {line.strip()}')

    # The same 3D page the HER2 run produced, with this target's numbers in it.
    page = None
    if hot:
        first, last = min(refres), max(refres)
        r4 = subprocess.run(
            [sys.executable, '-I', os.path.join(HERE, 'prepare_target.py'),
             os.path.join(cache, f'{ref_pdb}.cif'), ref_chain,
             '--window', f'{first}-{last}', '--hotspots', spec,
             '--uniprot', fasta, '--outdir', tdir, '--name', f'{stem}_full'],
            capture_output=True, text=True)
        whole = os.path.join(tdir, f'{stem}_full_target.pdb')
        if r4.returncode == 0 and os.path.exists(whole):
            import demo_page                                      # noqa: E402
            cdist = {k: v for k, v in A.get('counts_distinct', {}).items()}
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
                           ['residues measured', len(refres)],
                           ['numbering vs UniProt', f'{ref_off:+d}']],
                 'patch': {'label': f'residues {hot[0] + ref_off}-{hot[1] + ref_off} '
                                    f'(UniProt), {hot[2]} of {nd} binders',
                           'resi': list(range(hot[0], hot[1] + 1)),
                           'binders': hot[3][:8]}},
                whole, os.path.join(atlas_out, 'binders.pdb'))
            print(f'        3D page: {page}')

    print('        the measured data:')
    summary = {'target': args.name, 'uniprot': args.uniprot, 'full_name': fullname,
               'sequence_length': len(seq),
               'structures_listed': len(pdbs), 'structures_used': len(ok),
               'reference': f'{ref_pdb} chain {ref_chain}',
               'reference_offset': ref_off,
               'complexes': len(entries), 'units': nu, 'distinct_binders': nd,
               'numbering_conventions': conv, 'unit_formats': fmts,
               'domains_uniprot': [{'name': n, 'first': a, 'last': b}
                                   for n, a, b in doms],
               'most_targeted_patch': (
                   {'first_uniprot': hot[0] + ref_off, 'last_uniprot': hot[1] + ref_off,
                    'distinct_binders': hot[2], 'binders': hot[3]} if hot else None),
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
        json.dump({'schema': 'epitope-pipeline-snapshot/1',
                   'uniprot': {'accession': args.uniprot, 'full_name': fullname,
                               'sequence': seq, 'signal': sig,
                               'domains': [list(d) for d in doms],
                               'ecd': list(ecd) if ecd else None,
                               'pdbs': [{'id': p['id'],
                                         'range': list(p['range']) if p['range'] else None}
                                        for p in pdbs]},
                   'relatives': rel_seqs,
                   'entities': ent_doc},
                  fh, sort_keys=True, separators=(',', ':'))
        fh.write('\n')

    lock_out = os.path.join(out, 'inputs.lock.json')
    with open(lock_out, 'w') as fh:
        json.dump({'schema': 'epitope-pipeline-lock/2',
                   'generated': datetime.datetime.now(datetime.timezone.utc)
                                .replace(microsecond=0).isoformat(),
                   'target': args.name, 'uniprot': args.uniprot,
                   'sequence_sha256': seq_sha, 'sequence_length': len(seq),
                   'region': list(region) if region else None,
                   'limit': lock['limit'] if lock else args.limit,
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
          f'{nd} distinct binders out.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
