#!/bin/bash
#SBATCH -p sapphire
#SBATCH -J extract_keys_bc
#SBATCH --array=0-0%20
#SBATCH -c 1
#SBATCH --mem=6G
#SBATCH -t 08:00:00
#SBATCH -o %x.%A_%a.out
#SBATCH -e %x.%A_%a.err
set -euo pipefail
: "${RUN:?Set RUN (e.g. bc_6)}"
: "${WORK:?Set WORK (run workdir, e.g. /.../bc_6)}"
module load python
mamba activate barcode_amplicon_env_dev

CHUNKS=$WORK/chunks
PARTS=128

TASK=$(printf "%04d" ${SLURM_ARRAY_TASK_ID})
CHUNK_FASTQ=$CHUNKS/chunk_${TASK}.fastq
if [[ ! -s "$CHUNK_FASTQ" ]]; then
  echo "Missing chunk: $CHUNK_FASTQ"
  exit 1
fi

OUTDIR=$WORK/keys_chunks/chunk_${TASK}
mkdir -p "$OUTDIR"

SAMPLE_TSV=$OUTDIR/sample_list.tsv
printf "%s\t%s\n" "$RUN" "$CHUNK_FASTQ" > "$SAMPLE_TSV"

barcodeAmplicon-extract-keys "$SAMPLE_TSV" --out-dir "$OUTDIR" --parts "$PARTS"
