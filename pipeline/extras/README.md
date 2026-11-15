# Helper scripts from the HER2 analysis

These were written for the HER2 challenge. `run_pipeline.py` does not run them and the tests do not
cover them. Their thresholds are their own and differ from the pipeline's (a residue counts as
exposed at 0.25 here, at 0.10 in `surface_props.py`). Treat them as worked examples.

- `contacts.py`: which amino acids of the antigen touch a partner chain
- `numbering.py`: how a file's numbering differs from UniProt's (HER2 only)
- `sites.py`: rank candidate patches on one structure
- `paratope.py`: which amino acids of an antibody touch the antigen
- `format_fit.py`: the size and shape of a patch
- `pick_sites.py`: combine the evidence into site picks. It expects a trastuzumab footprint in the atlas.
