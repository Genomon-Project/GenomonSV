#! /bin/sh
#$ -S /bin/sh
#$ -cwd

INPUTDIR=$1
OUTPUTDIR=$2
REGION=$3

echo "samtools view ${INPUTDIR}/${REGION}.bam | perl bamPerseSV2.pl - > ${OUTPUTDIR}/${REGION}.improperPair.temp.txt"
samtools view ${INPUTDIR}/${REGION}.bam | perl bamPerseSV2.pl - > ${OUTPUTDIR}/${REGION}.improperPair.temp.txt

