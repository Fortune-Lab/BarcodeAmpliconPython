#!/bin/bash
set -euo pipefail

# Usage:
#   bash submit_pipeline.sh <RUN> <WORK> <FASTQ_R1_GZ> [PLOT|NONE]
#
# Example:
#   RUN=bc_6
#   WORK=/n/netscratch/.../Ehrt-VV-20974_2026_05_27/$RUN
#   FASTQ=/n/netscratch/.../bc_6_S375_L002_R1_001.fastq.gz
#   bash submit_pipeline.sh "$RUN" "$WORK" "$FASTQ" PLOT
#
# If 4th arg is "NONE", plot step is skipped.

RUN="${1:?RUN required (e.g. bc_6)}"
WORK="${2:?WORK required (e.g. /.../bc_6)}"
FASTQ="${3:?FASTQ required (full path to R1 fastq.gz)}"
DO_PLOT="${4:-PLOT}"

SPLIT_SBATCH="split_RUN.sbatch"
EXTRACT_SBATCH="extract_keys_RUN.sbatch"
MERGE_SBATCH="merge_keys_RUN.sbatch"
REDUCE_SBATCH="reduce_mcounts_RUN.sbatch"
COMBINE_SBATCH="combine_mcounts_RUN.sbatch"
THRESH_SBATCH="threshold_regression_RUN.sbatch"
PLOT_SBATCH="plot_threshold_RUN.sbatch"

CHUNKS_DIR="$WORK/chunks"

for s in "$SPLIT_SBATCH" "$EXTRACT_SBATCH" "$MERGE_SBATCH" "$REDUCE_SBATCH" "$COMBINE_SBATCH" "$THRESH_SBATCH"; do
  if [[ ! -f "$s" ]]; then
    echo "ERROR: missing sbatch script: $s" >&2
    exit 1
  fi
done
if [[ "$DO_PLOT" != "NONE" && ! -f "$PLOT_SBATCH" ]]; then
  echo "ERROR: plot requested but missing: $PLOT_SBATCH" >&2
  exit 1
fi

echo "RUN=$RUN"
echo "WORK=$WORK"
echo "FASTQ=$FASTQ"
echo "CHUNKS_DIR=$CHUNKS_DIR"
echo "DO_PLOT=$DO_PLOT"
echo

# 1) split (needs RUN/WORK/FASTQ)
jid_split=$(
  sbatch --parsable \
    --export=ALL,RUN="$RUN",WORK="$WORK",FASTQ="$FASTQ" \
    "$SPLIT_SBATCH"
)
echo "Submitted split: $jid_split ($SPLIT_SBATCH)"

echo
echo "Waiting for split to finish so we can count chunks..."
while true; do
  state=$(sacct -j "$jid_split" --format=State --noheader | head -n 1 | awk '{print $1}')
  if [[ "$state" == "COMPLETED" ]]; then
    echo "Split completed."
    break
  fi
  if [[ "$state" == "FAILED" || "$state" == "CANCELLED" || "$state" == "TIMEOUT" ]]; then
    echo "ERROR: Split job ended in state=$state" >&2
    exit 1
  fi
  sleep 20
done

# Count chunks (split created chunk_*.fastq)
N=$(ls -1 "$CHUNKS_DIR"/chunk_*.fastq 2>/dev/null | wc -l | awk '{print $1}')
if [[ "$N" -le 0 ]]; then
  echo "ERROR: No chunk_*.fastq found in $CHUNKS_DIR" >&2
  exit 1
fi
echo "Detected NCHUNKS=$N"

# 2) extract array (needs RUN/WORK; array bounds set here)
jid_extract=$(
  sbatch --parsable \
    --dependency=afterok:"$jid_split" \
    --array=0-$((N-1))%20 \
    --export=ALL,RUN="$RUN",WORK="$WORK" \
    "$EXTRACT_SBATCH"
)
echo "Submitted extract array: $jid_extract ($EXTRACT_SBATCH) [0..$((N-1))]"

# 3) merge
jid_merge=$(
  sbatch --parsable \
    --dependency=afterok:"$jid_extract" \
    --export=ALL,RUN="$RUN",WORK="$WORK" \
    "$MERGE_SBATCH"
)
echo "Submitted merge: $jid_merge ($MERGE_SBATCH)"

# 4) reduce (array over partitions is inside sbatch file)
jid_reduce=$(
  sbatch --parsable \
    --dependency=afterok:"$jid_merge" \
    --export=ALL,RUN="$RUN",WORK="$WORK" \
    "$REDUCE_SBATCH"
)
echo "Submitted reduce: $jid_reduce ($REDUCE_SBATCH)"

# 5) combine
jid_combine=$(
  sbatch --parsable \
    --dependency=afterok:"$jid_reduce" \
    --export=ALL,RUN="$RUN",WORK="$WORK" \
    "$COMBINE_SBATCH"
)
echo "Submitted combine: $jid_combine ($COMBINE_SBATCH)"

# 6) threshold
jid_thr=$(
  sbatch --parsable \
    --dependency=afterok:"$jid_combine" \
    --export=ALL,RUN="$RUN",WORK="$WORK" \
    "$THRESH_SBATCH"
)
echo "Submitted threshold: $jid_thr ($THRESH_SBATCH)"

# 7) plot (optional)
if [[ "$DO_PLOT" != "NONE" ]]; then
  jid_plot=$(
    sbatch --parsable \
      --dependency=afterok:"$jid_thr" \
      --export=ALL,RUN="$RUN",WORK="$WORK" \
      "$PLOT_SBATCH"
  )
  echo "Submitted plot: $jid_plot ($PLOT_SBATCH)"
else
  echo "Plot skipped."
fi

echo
echo "DONE submitting. Track with:"
echo "  squeue -u $USER"
echo "  sacct -j $jid_split,$jid_extract,$jid_merge,$jid_reduce,$jid_combine,$jid_thr --format=JobID,State,Elapsed"
