# How Micron Was Created

## Goal

Create a small, test-friendly DB subset (`micron/`) from full source DBs that preserves enough template/MSA content for ProteinFold test runs.

## Input Signals

- Archived per-sample hit artifacts from a known-good full-DB run (for this rebuild, `work_dev/Mar-16`).
- Sample list from CSV (for example `fulltest.csv`).
- Source DB root (for this rebuild, `/srv/scratch/sbf/dbs/proteinfold_dbs`).

## Core Strategy

1. Discover recent sample-specific hit files within a time window.
2. Extract top N hit IDs per sample (default N=4).
3. Aggregate IDs across samples and cap total set size (for example ~12).
4. Rewrite selected micron DB files from source using those IDs.
5. Reconcile template indices (`pdb70`, `pdb100`) with mmCIF files.
6. Sync model parameter assets into `params/`, including Boltz checkpoints and chemistry files.
7. Prune unreferenced mmCIF files to `pdb_mmcif/orphaned/`.

## Updated Paths in `micron/`

- `uniref90/uniref90.fasta`
- `mgnify/mgy_clusters.fa`
- `small_bfd/bfd-first_non_consensus_sequences.fasta`
- `pdb70/pdb70_{a3m,cs219,hhm}.{ffdata,ffindex}`
- `pdb100/pdb100_2021Mar03_{a3m,cs219,hhm,pdb}.{ffdata,ffindex}`
- `pdb_mmcif/mmcif_files/*.cif`
- `params/`
- `params/mols/`

## Safety Rules

- If ID filtering yields empty FASTA output, fallback to top N source records.
- Required mmCIFs are derived from truncated template indices.
- Missing required mmCIFs trigger a hard error.
- Extra mmCIFs are moved to `pdb_mmcif/orphaned/` for review.
- `params/mols/` must retain the canonical Boltz molecule `.pkl` set required for inference. Additional `.pkl` files are optional.
