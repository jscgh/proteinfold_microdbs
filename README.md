# Micron Reproducibility Repo

This repository captures how the `micron/` database subset is produced from ProteinFold DBs and recent MSA/template hits.

## Contents

- `scripts/truncate_micron_from_sample_hits.sh`
  - Main reproducible truncation script.
  - Uses per-sample recent hit artifacts, aggregates top hits, updates DB subsets, reconciles `pdb70`/`pdb100` with `pdb_mmcif/mmcif_files`.
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
- It logs run output and prints a change summary of updated files.
