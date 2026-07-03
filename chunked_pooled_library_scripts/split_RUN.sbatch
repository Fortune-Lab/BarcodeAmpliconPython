#!/bin/bash
#SBATCH -p sapphire
#SBATCH -J split_bc
#SBATCH -c 8
#SBATCH --mem=8G
#SBATCH -t 02:00:00
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
set -euo pipefail
: "${RUN:?Set RUN (e.g. bc_6)}"
: "${WORK:?Set WORK (run workdir, e.g. /.../bc_6)}"
: "${FASTQ:?Set FASTQ (full path to R1 fastq.gz)}"
CHUNKS=$WORK/chunks
READS_PER_CHUNK=10000000

mkdir -p "$CHUNKS" "$WORK/logs"
cd "$CHUNKS"

pigz -dc -p "${SLURM_CPUS_PER_TASK:-1}" "$FASTQ" \
  | split -l $((4*READS_PER_CHUNK)) -d -a 4 --additional-suffix=.fastq - chunk_

echo "NCHUNKS=$(ls -1 chunk_*.fastq | wc -l)"
