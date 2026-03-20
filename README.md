# Micron Reproducibility Repo

This repository captures how the `micron/` database subset is produced from ProteinFold full DBs and archived MSA/template hits.

## Contents

- `scripts/truncate_micron_from_sample_hits.sh`
  - Main reproducible truncation script.
  - Uses per-sample hit artifacts from archived run directories, aggregates top hits, updates DB subsets, reconciles `pdb70`/`pdb100` with `pdb_mmcif/mmcif_files`, and syncs `params/` for Boltz.
- `scripts/submit_truncate_job.pbs`
  - PBS batch wrapper for running the main truncation script on compute.
- `docs/micron_creation.md`
  - Design and provenance notes.
- `docs/repro_steps.md`
  - Step-by-step reproducibility procedure.

## Quick Start (Compute Job)

Edit `scripts/submit_truncate_job.pbs` variables, then submit:

```bash
qsub scripts/submit_truncate_job.pbs
```

## Notes

- The main script is designed to avoid producing unusable template DB outputs.
- For the March 16 rebuild, use `/srv/scratch/z5378336/inputs/test_inputs/work_dev/Mar-16` with `WINDOW_MIN=0` so archived files are included.
- The script now populates `micron/params` from the full params tree and always keeps the canonical Boltz `params/mols/*.pkl` set required for inference. `BOLTZ_MOLS_KEEP` adds extra molecule `.pkl` files on top of that baseline.
- It logs run output and prints a change summary of updated files.
