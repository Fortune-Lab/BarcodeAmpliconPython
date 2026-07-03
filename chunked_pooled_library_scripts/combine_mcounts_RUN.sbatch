#!/bin/bash
#SBATCH -p sapphire
#SBATCH -J combine_mcounts_bc
#SBATCH -c 1
#SBATCH --mem=2G
#SBATCH -t 00:30:00
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
set -euo pipefail
: "${RUN:?Set RUN (e.g. bc_6)}"
: "${WORK:?Set WORK (run workdir, e.g. /.../bc_6)}"
OUT=$WORK/${RUN}_mcounts.tsv

echo -e "run\tqtag_id\tbarcode\tmcount" > "$OUT"
for f in $(ls -1 $WORK/${RUN}_mcounts_part_*.tsv | sort); do
  tail -n +2 "$f" >> "$OUT"
done

wc -l "$OUT"
head -n 5 "$OUT"
