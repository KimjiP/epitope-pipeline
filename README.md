# Epitope pipeline

**Where should a new antibody bind?** Give this one UniProt accession and it answers from
every published structure of that antigen: who already binds it, where, and which surface
the field keeps returning to.

```bash
python3 -I pipeline/run_pipeline.py --uniprot P00533 --name EGFR --out runs/egfr
```

About ten seconds. Python standard library only — no install, no dependencies, no account.

---

## The six steps

The whole method is one pass. Each step is a script you can also run on its own.

| | Step | What it does | Script |
|---|---|---|---|
| **1** | Ask what is already known | Read the antigen's sequence, its domain boundaries and its list of structures from the UniProt entry | `run_pipeline.py` |
| **2** | Collect every published structure | Keep the ones covering the part a drug can reach, and download them from the PDB | `run_pipeline.py` |
| **3** | Put them all in one coordinate frame | Recover each file's own residue numbering by alignment, then superpose every complex onto one reference | `atlas.py` |
| **4** | Count distinct binders, not files | One antibody solved six times is one binder. Split complexes into binding units, merge units by sequence identity | `atlas.py` |
| **5** | Measure the surface and score every site | Exposure, hydrophobicity, burial and contact density per residue; which properties predict where binders land | `surface_props.py`, `pick_sites.py` |
| **6** | Hand over what a design model needs | A hotspot list verified against the structure, a clean target file, and a 3D page | `prepare_target.py`, `demo_page.py` |

### Why steps 3 and 4 are the hard ones

**Nothing in a PDB file says which residue numbering it uses.** Across one antigen's
structures you will find the mature numbering, the UniProt numbering and construct offsets,
with no field that distinguishes them. Step 3 recovers each file's offset by aligning its
chains to the reference sequence. Getting this wrong moves every epitope by a constant and
produces a plausible, wrong answer.

**Counting files is not counting binders.** A popular antibody appears in many entries, and
one entry can hold two different antibodies. Step 4 splits each complex into binding units
(chains that touch each other), then merges units whose variable domains are at least 95%
identical. On HER2 this is the difference between 36 units and 24 actual binders, and it
reverses which epitope looks most popular.

---

## Reproducibility

Every run writes two files that pin it:

- **`inputs.lock.json`**: the exact PDB entries used, the SHA-256 of each file as downloaded,
  the reference structure chosen, the relatives excluded, and a hash of the UniProt sequence
  everything was numbered against.
- **`inputs.snapshot.json`**: everything else the run reads from the network. That is the
  UniProt facts (sequence, domains, topology), the sequences of the excluded relatives, and the
  PDB's chain names and resolutions. The lock holds its hash.

```bash
python3 -I pipeline/run_pipeline.py --from-lock examples/egfr/inputs.lock.json --out runs/egfr-replay
```

A replay uses that exact list, that exact reference and that exact metadata, and checks every
file against its recorded hash. Once the structure files are cached it needs **no network at
all**. If a structure has been re-deposited it says so, instead of quietly producing a
different answer.

Without the lock a rerun is not the same run: the PDB grows, so `--limit` takes a different
set, and the reference the scan picks can change with it.

A run only ever reads the entries it downloaded, named explicitly, so one `--cache` shared
between antigens cannot leak another antigen's structures into the atlas.

**Check the pipeline still reproduces its examples:**

```bash
./reproduce.sh
```

This replays both locks, then replays them again with the network blocked, and compares each
result against the committed `summary.json`. It prints a line per number and exits non-zero on
any difference. It downloads about 110 MB the first time and caches it.

Replays on Python 3.9, 3.13 and 3.14 give identical result files. Only timestamps and the run
folder's own path differ.

**Check the maths without a network:**

```bash
python3 -I pipeline/atlas.py --selftest
```

Superposing a point set onto itself must give the identity matrix and exactly zero RMSD; a
known rotation must come back exactly; the aligner must score 1.000 for identical sequences
and below 0.5 for unrelated ones. This runs on every push.

---

## Options

| Option | |
|---|---|
| `--uniprot ACC` | the antigen. Required, unless `--from-lock` supplies it. |
| `--name NAME` | a short name for the output files. Required, unless `--from-lock` supplies it. |
| `--limit N` | how many structures to download (default 30). The first N of those that pass the region filter. |
| `--region LO-HI` | keep only structures covering these residues. Default: the extracellular part, from the UniProt entry's own topology annotation. |
| `--exclude-uniprot A,B` | close relatives of the antigen. Their chains bind as receptor partners, not as drugs, and counting them as binders distorts the result. |
| `--from-lock PATH` | replay a previous run exactly, offline once the structures are cached. |
| `--out DIR` | where everything goes. |
| `--cache DIR` | where structure files live (default `OUT/pdb-cache`). Share one across runs. |

### What a run writes

| File | |
|---|---|
| `inputs.lock.json` | the exact inputs, for replay |
| `inputs.snapshot.json` | the network data the run read, so a replay needs no network |
| `summary.json` | every headline number |
| `atlas-out/atlas.json` | each binder, its epitope, its cluster, per-residue counts |
| `atlas-out/binders.pdb` | every binder superposed into one frame |
| `props.json` | per-residue exposure, hydrophobicity, burial |
| `target/<name>_hotspots.txt` | the hotspot list, every residue letter verified against the structure |
| `viewer/index.html` | a 3D page, the surface coloured by binder density. Serve it with `python3 -m http.server --directory runs/egfr/viewer`; a browser will not load its data files from a `file://` page |
| `atlas.log` | the full atlas run, including everything skipped and why |

---

## Examples

`examples/her2` and `examples/egfr` hold a lock, a snapshot, a summary and a saved log for each. They are the
pipeline's output, not its purpose — the point is that you can run it on your own antigen and
on whatever is in the PDB today.

| | HER2 | EGFR |
|---|---|---|
| Structures on the UniProt entry | 63 | 388 |
| Used — something bound | 31 | 33 |
| Binding units | 36 | 37 |
| Distinct binders | 24 | 19 |
| Most-targeted patch (UniProt numbering) | 167–190 | 435–458 |
| Reference the pipeline chose | 8PWH chain E | 3NJP chain A |
| That file's numbering, against UniProt | +22 | +24 |

---

## Limits

- **It reads the first N structures that pass the region filter, not the best N.** Raise
  `--limit` for a thorough run.
- **Not everything that binds is a drug.** On EGFR the pipeline finds the antigen's own growth
  factors; on PD-L1 it finds PD-1. It reports them as not antibody-shaped rather than counting
  them as competitors.
- **Binder names and formats come from the depositors' text**, so read the names, not the
  labels.
- **The most-targeted patch is a starting point for a person**, not a site choice: a
  24-residue window ranked by how many distinct binders reach into it.

---

## Requirements

Python 3.9 or newer. Nothing else. Network access for UniProt (`rest.uniprot.org`) and the
PDB (`files.rcsb.org`, `data.rcsb.org`) while a run downloads; a replay from cache needs none.

## Data

Structures from the [RCSB PDB](https://www.rcsb.org/). Sequences, domain and topology
annotations from [UniProt](https://www.uniprot.org/) (CC BY 4.0). Hydrophobicity scale from
Kyte and Doolittle 1982; maximum accessible surface areas from Tien et al. 2013; superposition
by Horn's quaternion method.

## Licence

MIT. See [LICENSE](LICENSE).
