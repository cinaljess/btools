#!/bin/bash
RUN_DIR=$1
SAMPLE=$2
SAMPLE2=$3
PREFIX=${SAMPLE%.fastq}
REFERENCE=$4
VARGS="DP,DP4,ADF,ADR,SP"
if [ -z $5 ]
  then
    source /broad/software/scripts/useuse
    reuse Python-2.7
    reuse Samtools
    reuse .bwa-0.7.10
    reuse .bedtools-2.17.0
    reuse R-3.0
    reuse Java-1.8
    reuse GATK3
    reuse Picard-Tools
    #smalt map -n 4 $RUN_DIR/$REFERENCE $RUN_DIR/$SAMPLE > $RUN_DIR/$PREFIX.sam

    samtools view -hbF -4 $RUN_DIR/$PREFIX.sam > $RUN_DIR/$PREFIX.mapped.bam

    samtools sort -O bam -o $RUN_DIR/$PREFIX.mapped.sorted.bam -T $RUN_DIR/$PREFIX.ttt $RUN_DIR/$PREFIX.mapped.bam

    samtools index -b $RUN_DIR/$PREFIX.mapped.sorted.bam

    samtools mpileup -u -d 1000 -L 1000 --VCF --output-tags ${VARGS} -f $RUN_DIR/$REFERENCE.fa -o $RUN_DIR/$PREFIX.mapped.indels $RUN_DIR/$PREFIX.mapped.sorted.bam

    samtools mpileup -u --skip-indels -d 1000 --VCF --output-tags ${VARGS} -f $RUN_DIR/$REFERENCE.fa -o $RUN_DIR/$PREFIX.mapped.mutations $RUN_DIR/$PREFIX.mapped.sorted.bam
  else
    echo "$RUN_DIR - $SAMPLE - $SAMPLE2 - $PREFIX - $REFERENCE"
    echo "smalt map -n 4 $RUN_DIR/$REFERENCE $RUN_DIR/$SAMPLE > $RUN_DIR/$PREFIX.sam"
    echo "samtools view -hbF -4 $RUN_DIR/$PREFIX.sam > $RUN_DIR/$PREFIX.mapped.bam"
    echo "samtools sort -O bam -o $RUN_DIR/$PREFIX.mapped.sorted.bam -T $RUN_DIR/$PREFIX.ttt $RUN_DIR/$PREFIX.mapped.bam"
    echo "samtools index -b $RUN_DIR/$PREFIX.mapped.sorted.bam"
    echo "samtools mpileup -u -d 1000 -L 1000 --VCF --output-tags $VARGS -f $RUN_DIR/$REFERENCE.fa -o $RUN_DIR/$PREFIX.mapped.indels $RUN_DIR/$PREFIX.mapped.sorted.bam"
    echo "samtools mpileup -u --skip-indels -d 1000 --VCF --output-tags $VARGS -f $RUN_DIR/$REFERENCE.fa -o $RUN_DIR/$PREFIX.mapped.mutations $RUN_DIR/$PREFIX.mapped.sorted.bam"
fi
