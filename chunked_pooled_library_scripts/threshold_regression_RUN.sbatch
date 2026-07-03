#!/bin/bash
#SBATCH -p sapphire
#SBATCH -J thr_reg_bc
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 02:00:00
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
set -euo pipefail
: "${RUN:?Set RUN (e.g. bc_6)}"
: "${WORK:?Set WORK (run workdir, e.g. /.../bc_6)}"
module load python
mamba activate barcode_amplicon_env_dev


INFILE=$WORK/${RUN}_mcounts.tsv
WINDOW=500
OUTFILE=$WORK/${RUN}_threshold_regression_window_${WINDOW}.tsv

barcodeAmplicon-threshold-regression \
  --in "$INFILE" \
  --window "$WINDOW" \
  --out "$OUTFILE" \
  --run-mode FALSE \
  --percent-cutoff 0
