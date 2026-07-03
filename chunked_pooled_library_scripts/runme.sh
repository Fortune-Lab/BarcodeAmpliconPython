RUN=bc_6
WORK=/n/netscratch/sfortune_lab/Lab/mchase/DUC_Cornell/Ehrt-VV-20974_2026_05_27/chunked_pooled_library_scripts/$RUN
FASTQ=/n/netscratch/sfortune_lab/Lab/mchase/DUC_Cornell/Ehrt-VV-20974_2026_05_27/bc_6_S375_L002_R1_001.fastq.gz

bash submit_pipeline.sh "$RUN" "$WORK" "$FASTQ" PLOT
