#!/bin/bash
#SBATCH -p sapphire
#SBATCH -J merge_keys_bc
#SBATCH -c 2
#SBATCH --mem=4G
#SBATCH -t 02:00:00
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
set -euo pipefail
: "${RUN:?Set RUN (e.g. bc_6)}"
: "${WORK:?Set WORK (run workdir, e.g. /.../bc_6)}"
PARTS=128
CHUNK_KEYS_BASE=$WORK/keys_chunks
MERGED=$WORK/keys_merged

mkdir -p "$MERGED"

# fail early if no chunk outputs
if ! ls -d "$CHUNK_KEYS_BASE"/chunk_* >/dev/null 2>&1; then
  echo "No chunk dirs under $CHUNK_KEYS_BASE"
  exit 1
fi

for p in $(seq 0 $((PARTS-1))); do
  P=$(printf "%03d" $p)
  out="$MERGED/keys_part_${P}.tsv"
  : > "$out"
  for d in "$CHUNK_KEYS_BASE"/chunk_*; do
    f="$d/$RUN/keys_part_${P}.tsv"
    [[ -s "$f" ]] && cat "$f" >> "$out"
  done
done

ls -lh "$MERGED"/keys_part_000.tsv "$MERGED"/keys_part_127.tsv
