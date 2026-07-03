#!/bin/bash
#SBATCH -p sapphire
#SBATCH -J reduce_mcounts_bc
#SBATCH --array=0-127%16
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 08:00:00
#SBATCH -o %x.%A_%a.out
#SBATCH -e %x.%A_%a.err
set -euo pipefail
: "${RUN:?Set RUN (e.g. bc_6)}"
: "${WORK:?Set WORK (run workdir, e.g. /.../bc_6)}"
module load python
mamba activate barcode_amplicon_env_dev

KEYS_DIR=$WORK/keys_merged

PART=${SLURM_ARRAY_TASK_ID}
PART_TAG=$(printf "%03d" $PART)
OUT=$WORK/${RUN}_mcounts_part_${PART_TAG}.tsv

barcodeAmplicon-reduce-mcounts \
  --keys-dir "$KEYS_DIR" \
  --run "$RUN" \
  --part "$PART" \
  --out "$OUT" \
  --tmp-dir "${SLURM_TMPDIR:-/tmp}" \
  --sort-mem 4G
