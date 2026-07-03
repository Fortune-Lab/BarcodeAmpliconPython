#!/bin/bash
#SBATCH -p sapphire
#SBATCH -J plot_bc
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 00:30:00
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
set -euo pipefail
: "${RUN:?Set RUN (e.g. bc_6)}"
: "${WORK:?Set WORK (run workdir, e.g. /.../bc_6)}"
module load python
mamba activate barcode_amplicon_env_dev

INFILE=$WORK/${RUN}_threshold_regression_window_500.tsv

barcodeAmplicon-plot-threshold \
  --outdir "$WORK" \
  --stride 2 \
  --xmax 50000 \
  --left_axis log10_counts \
  "$INFILE"
