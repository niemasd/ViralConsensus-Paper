#!/usr/bin/env bash
# subsample FASTQ pair w/ sektk, map reads w/ Minimap2, and run both iVar and AmpliPy pipelines
if [ "$#" -ne 7 ] ; then
    echo "USAGE: $0 <FASTQ> <TOT_NUM_READS> <REF_GENOME_FAS> <REF_GENOME_MMI> <PRIMER_BED> <OUT_ZIP> <RNG_SEED>"; exit 1
fi

# parse and check args
FQ=$1 ; N=$2 ; REF_FAS=$3 ; REF_MMI=$4 ; PRIMER_BED=$5 ; OUT_ZIP=$6 ; SEED=$7 ; TMP_OUT_DIR=$(mktemp -d)
if [ ! -f "$FQ" ] ; then
    echo "File not found: $FQ" ; exit 1
elif [ -f "$OUT_ZIP" ] ; then
    echo "File already exists: $OUT_ZIP" ; exit 1
fi

# subsample FASTQ pair using seqtk
FQ_SUB="$TMP_OUT_DIR/sub.fastq.gz"
TIME_SEQTK="$TMP_OUT_DIR/time.01.seqtk.subsample.txt"
/usr/bin/time -v -o "$TIME_SEQTK" seqtk sample "-s$SEED" "$FQ" $N | gzip -9 > "$FQ_SUB"

# map reads using Minimap2
UNTRIMMED_BAM="$TMP_OUT_DIR/untrimmed.bam"
TIME_MINIMAP2="$TMP_OUT_DIR/time.02.minimap2.txt"
LOG_MINIMAP2="$TMP_OUT_DIR/log.minimap2.txt"
/usr/bin/time -v -o "$TIME_MINIMAP2" minimap2 -a -x sr "$REF_MMI" "$FQ_SUB" 2> "$LOG_MINIMAP2" | samtools view -S -b > "$UNTRIMMED_BAM"

# sort untrimmed BAM using samtools
UNTRIMMED_SORTED_BAM="$TMP_OUT_DIR/untrimmed.sorted.bam"
TIME_SAMTOOLS_SORT_UNTRIMMED="$TMP_OUT_DIR/time.03.samtools.sort.untrimmed.txt"
/usr/bin/time -v -o "$TIME_SAMTOOLS_SORT_UNTRIMMED" samtools sort -o "$UNTRIMMED_SORTED_BAM" "$UNTRIMMED_BAM"

# trim reads using iVar Trim
IVAR_TRIMMED_BAM_PREFIX="$TMP_OUT_DIR/trimmed.ivar"
IVAR_TRIMMED_BAM="$IVAR_TRIMMED_BAM_PREFIX.bam"
TIME_IVAR_TRIM="$TMP_OUT_DIR/time.04.ivar.trim.txt"
LOG_IVAR_TRIM="$TMP_OUT_DIR/log.ivar.trim.txt"
/usr/bin/time -v -o "$TIME_IVAR_TRIM" ivar trim -x 5 -e -i "$UNTRIMMED_SORTED_BAM" -b "$PRIMER_BED" -p "$IVAR_TRIMMED_BAM_PREFIX" > "$LOG_IVAR_TRIM" 2>&1

# sort iVar-trimmed BAM using samtools
IVAR_TRIMMED_SORTED_BAM="$TMP_OUT_DIR/trimmed.ivar.sorted.bam"
TIME_SAMTOOLS_SORT_IVAR_TRIMMED="$TMP_OUT_DIR/time.05.samtools.sort.ivar.trimmed.txt"
/usr/bin/time -v -o "$TIME_SAMTOOLS_SORT_IVAR_TRIMMED" samtools sort -o "$IVAR_TRIMMED_SORTED_BAM" "$IVAR_TRIMMED_BAM"

# generate pile-up using samtools
PILEUP="$TMP_OUT_DIR/trimmed.ivar.sorted.pileup.txt"
TIME_SAMTOOLS_PILEUP="$TMP_OUT_DIR/time.06.samtools.pileup.txt"
LOG_SAMTOOLS_PILEUP="$TMP_OUT_DIR/log.samtools.pileup.txt"
/usr/bin/time -v -o "$TIME_SAMTOOLS_PILEUP" samtools mpileup -A -aa -d 0 -Q 0 --reference "$REF_FAS" "$IVAR_TRIMMED_SORTED_BAM" > "$PILEUP" 2> "$LOG_SAMTOOLS_PILEUP"

# call variants using iVar Variants
IVAR_VARIANTS_TSV="$TMP_OUT_DIR/trimmed.ivar.sorted.pileup.variants.ivar.tsv"
TIME_IVAR_VARIANTS="$TMP_OUT_DIR/time.07.ivar.variants.txt"
LOG_IVAR_VARIANTS="$TMP_OUT_DIR/log.ivar.variants.txt"
/usr/bin/time -v -o "$TIME_IVAR_VARIANTS" bash -c "cat $PILEUP | ivar variants -r $REF_FAS -p $IVAR_VARIANTS_TSV -m 10 > $LOG_IVAR_VARIANTS 2>&1"

# call consensus using iVar Consensus
IVAR_CONSENSUS_PREFIX="$TMP_OUT_DIR/trimmed.ivar.sorted.pileup.consensus.ivar"
TIME_IVAR_CONSENSUS="$TMP_OUT_DIR/time.08.ivar.consensus.txt"
LOG_IVAR_CONSENSUS="$TMP_OUT_DIR/log.ivar.consensus.txt"
/usr/bin/time -v -o "$TIME_IVAR_CONSENSUS" bash -c "cat $PILEUP | ivar consensus -p $IVAR_CONSENSUS_PREFIX -m 10 -n N -t 0.5 > $LOG_IVAR_CONSENSUS 2>&1"

# trim reads using AmpliPy Trim
AMPLIPY_TRIMMED_BAM="$TMP_OUT_DIR/trimmed.amplipy.bam"
TIME_AMPLIPY_TRIM="$TMP_OUT_DIR/time.09.amplipy.trim.txt"
LOG_AMPLIPY_TRIM="$TMP_OUT_DIR/log.amplipy.trim.txt"
/usr/bin/time -v -o "$TIME_AMPLIPY_TRIM" AmpliPy.py trim -i "$UNTRIMMED_SORTED_BAM" -p "$PRIMER_BED" -r "$REF_FAS" -o "$AMPLIPY_TRIMMED_BAM" -x 5 -e 2> "$LOG_AMPLIPY_TRIM"

# call variants using AmpliPy Variants
AMPLIPY_VARIANTS_VCF="$TMP_OUT_DIR/trimmed.ivar.sorted.variants.amplipy.vcf"
TIME_AMPLIPY_VARIANTS="$TMP_OUT_DIR/time.10.amplipy.variants.txt"
LOG_AMPLIPY_VARIANTS="$TMP_OUT_DIR/log.amplipy.variants.txt"
/usr/bin/time -v -o "$TIME_AMPLIPY_VARIANTS" AmpliPy.py variants -i "$IVAR_TRIMMED_SORTED_BAM" -r "$REF_FAS" -o "$AMPLIPY_VARIANTS_VCF" -md 10 2> "$LOG_AMPLIPY_VARIANTS"

# call consensus using AmpliPy Consensus
AMPLIPY_CONSENSUS_FAS="$TMP_OUT_DIR/trimmed.ivar.sorted.consensus.amplipy.fas"
TIME_AMPLIPY_CONSENSUS="$TMP_OUT_DIR/time.11.amplipy.consensus.txt"
LOG_AMPLIPY_CONSENSUS="$TMP_OUT_DIR/log.amplipy.consensus.txt"
/usr/bin/time -v -o "$TIME_AMPLIPY_CONSENSUS" AmpliPy.py consensus -i "$IVAR_TRIMMED_SORTED_BAM" -r "$REF_FAS" -o "$AMPLIPY_CONSENSUS_FAS" -md 10 -n N -mf 0.5 2> "$LOG_AMPLIPY_CONSENSUS"

# trim reads + call variants + call consensus using AmpliPy AIO
AMPLIPY_AIO_TRIMMED_BAM="$TMP_OUT_DIR/trimmed.amplipy.aio.bam"
AMPLIPY_AIO_VARIANTS_VCF="$TMP_OUT_DIR/trimmed.amplipy.aio.variants.vcf"
AMPLIPY_AIO_CONSENSUS_FAS="$TMP_OUT_DIR/trimmed.amplipy.aio.consensus.fas"
TIME_AMPLIPY_AIO="$TMP_OUT_DIR/time.12.amplipy.aio.txt"
LOG_AMPLIPY_AIO="$TMP_OUT_DIR/log.amplipy.aio.txt"
/usr/bin/time -v -o "$TIME_AMPLIPY_AIO" AmpliPy.py aio -i "$UNTRIMMED_SORTED_BAM" -p "$PRIMER_BED" -r "$REF_FAS" -ot "$AMPLIPY_AIO_TRIMMED_BAM" -ov "$AMPLIPY_AIO_VARIANTS_VCF" -oc "$AMPLIPY_AIO_CONSENSUS_FAS" -x 5 -e -mdv 10 -mdc 10 -n N -mfc 0.5 2> "$LOG_AMPLIPY_AIO"

# call consensus using ViralConsensus
VIRALCONSENSUS_CONSENSUS_FAS="$TMP_OUT_DIR/consensus.viralconsensus.fas"
TIME_VIRALCONSENSUS_CONSENSUS="$TMP_OUT_DIR/time.13.viralconsensus.txt"
LOG_AMPLIPY_CONSENSUS="$TMP_OUT_DIR/log.viralconsensus.txt"
/usr/bin/time -v -o "$TIME_VIRALCONSENSUS_CONSENSUS" viral_consensus -i "$UNTRIMMED_BAM" -p "$PRIMER_BED" -r "$REF_FAS" -o "$VIRALCONSENSUS_CONSENSUS_FAS" -op /dev/null -oi /dev/null -po 5

# zip output and clean up
rm -f $TMP_OUT_DIR/*.fastq.gz $TMP_OUT_DIR/*.bam $TMP_OUT_DIR/*.bam.bai $TMP_OUT_DIR/*.pileup.txt  # delete FASTQs and BAMs to save space
zip -j -9 "$OUT_ZIP" $TMP_OUT_DIR/*
rm -rf "$TMP_OUT_DIR"
