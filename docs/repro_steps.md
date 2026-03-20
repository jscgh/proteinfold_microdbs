# Reproducibility Steps

## 1. Preconditions

- Full source DB tree available (for example `/srv/scratch/sbf/dbs/proteinfold_dbs`).
- Work directory with archived MSA/template artifacts.
- Sample CSV available.
- Compute environment can load required modules.

## 2. Recommended Run (from compute)

```bash
scripts/truncate_micron_from_sample_hits.sh \
  /srv/scratch/z5378336/inputs/test_inputs/work_dev/Mar-16 \
  /srv/scratch/sbf/dbs/proteinfold_dbs \
  /srv/scratch/sbf/dbs/proteinfold_microdbs/micron \
  /srv/scratch/z5378336/inputs/test_inputs/fulltest.csv \
  4 \
  12 \
  /srv/scratch/sbf/dbs/proteinfold_microdbs/micron/truncate_micron_from_sample_hits.log \
  0 \
  12 \
  /srv/scratch/sbf/dbs/proteinfold_dbs/params \
  0
```

Argument order:

1. `WORK_DIR`
2. `SOURCE_DB_DIR`
3. `MICRON_DIR`
4. `SAMPLES_CSV`
5. `N_PER_SAMPLE`
6. `MAX_TOTAL`
7. `LOG_FILE`
8. `WINDOW_MIN`
9. `COLABFOLD_N`
10. `PARAMS_SOURCE_DIR`
11. `BOLTZ_MOLS_KEEP`

## 3. PBS Submission

Use `scripts/submit_truncate_job.pbs` and submit with `qsub`.

## 4. Validation Checklist

- Log contains `Change summary` entries.
- `micron/pdb70/pdb70_a3m.ffindex` and `micron/pdb100/pdb100_2021Mar03_a3m.ffindex` are non-empty.
- For every required template ID in those indices, matching `micron/pdb_mmcif/mmcif_files/<pdbid>.cif` exists.
- FASTA outputs are non-empty.
- `micron/params/boltz1_conf.ckpt`, `micron/params/boltz2_aff.ckpt`, and `micron/params/boltz2_conf.ckpt` exist.
- `micron/params/mols/` contains the canonical Boltz molecule `.pkl` set required for inference, plus any extra `.pkl` files requested via `BOLTZ_MOLS_KEEP`.
