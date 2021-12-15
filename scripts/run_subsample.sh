#!/usr/bin/env bash
# subsample FASTQ pair w/ sektk, map reads w/ Minimap2, and run both iVar and AmpliPy pipelines
if [ "$#" -ne 7 ] ; then
    echo "USAGE: $0 <R1_FASTQ> <R2_FASTQ> <TOT_NUM_READS> <REF_GENOME_MMI> <PRIMER_BED> <OUT_ZIP> <RNG_SEED>"; exit 1
fi

# parse and check args
R1_FQ=$1 ; R2_FQ=$2 ; N=$3 ; REF_MMI=$4 ; PRIMER_BED=$5 ; OUT_ZIP=$6 ; SEED=$7 ; N_OVER_2=$(($N/2)); TMP_OUT_DIR=$(mktemp -d)
if [ ! -f "$R1_FQ" ] ; then
    echo "File not found: $R1_FQ" ; exit 1
elif [ ! -f "$R2_FQ" ] ; then
    echo "File not found: $R2_FQ" ; exit 1
elif [ -f "$OUT_ZIP" ] ; then
    echo "File already exists: $OUT_ZIP" ; exit 1
fi

# subsample FASTQ pair using seqtk
R1_FQ_SUB="$TMP_OUT_DIR/sub.R1.fastq.gz"
R2_FQ_SUB="$TMP_OUT_DIR/sub.R2.fastq.gz"
TIME_SEQTK_R1="$TMP_OUT_DIR/time.01.seqtk.subsample.R1.txt"
TIME_SEQTK_R2="$TMP_OUT_DIR/time.02.seqtk.subsample.R2.txt"
/usr/bin/time -v -o "$TIME_SEQTK_R1" seqtk sample "-s$SEED" "$R1_FQ" $N_OVER_2 | gzip -9 > "$R1_FQ_SUB"
/usr/bin/time -v -o "$TIME_SEQTK_R2" seqtk sample "-s$SEED" "$R2_FQ" $N_OVER_2 | gzip -9 > "$R2_FQ_SUB"

# map reads using Minimap2
UNTRIMMED_BAM="$TMP_OUT_DIR/untrimmed.bam"
TIME_MINIMAP2="$TMP_OUT_DIR/time.03.minimap2.txt"
LOG_MINIMAP2="$TMP_OUT_DIR/log.minimap2.txt"
/usr/bin/time -v -o "$TIME_MINIMAP2" minimap2 -a -x sr "$REF_MMI" "$R1_FQ_SUB" "$R2_FQ_SUB" 2> "$LOG_MINIMAP2" | samtools view -S -b > "$UNTRIMMED_BAM"
#rm -f "$R1_FQ_SUB" "$R2_FQ_SUB" # delete subsampled FQs (do this if they get too big)

# sort untrimmed BAM using samtools
UNTRIMMED_SORTED_BAM="$TMP_OUT_DIR/untrimmed.sorted.bam"
TIME_SAMTOOLS_SORT_UNTRIMMED="$TMP_OUT_DIR/time.04.samtools.sort.untrimmed.txt"
/usr/bin/time -v -o "$TIME_SAMTOOLS_SORT_UNTRIMMED" samtools sort -o "$UNTRIMMED_SORTED_BAM"

# zip output and clean up
zip -j -9 "$OUT_ZIP" $TMP_OUT_DIR/*
rm -rf "$TMP_OUT_DIR"
