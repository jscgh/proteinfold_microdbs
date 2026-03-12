# Reproducibility Steps

## 1. Preconditions

- Full source DB tree available (for example `old/`).
- Work directory with recent MSA/template artifacts.
- Sample CSV available.
- Compute environment can load required modules.

## 2. Recommended Run (from compute)

```bash
scripts/truncate_micron_from_sample_hits.sh \
  /srv/scratch/z5378336/inputs/test_inputs/work/Mar-11 \
  /srv/scratch/sbf/dbs/proteinfold_microdbs/old \
  /srv/scratch/sbf/dbs/proteinfold_microdbs/micron \
  /srv/scratch/z5378336/inputs/test_inputs/fulltest.csv \
  4 \
  12 \
  /srv/scratch/sbf/dbs/proteinfold_microdbs/micron/truncate_micron_from_sample_hits.log \
  70
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

## 3. PBS Submission

Use `scripts/submit_truncate_job.pbs` and submit with `qsub`.

## 4. Validation Checklist

- Log contains `Change summary` entries.
- `micron/pdb70/pdb70_a3m.ffindex` and `micron/pdb100/pdb100_2021Mar03_a3m.ffindex` are non-empty.
- For every required template ID in those indices, matching `micron/pdb_mmcif/mmcif_files/<pdbid>.cif` exists.
- FASTA outputs are non-empty.
