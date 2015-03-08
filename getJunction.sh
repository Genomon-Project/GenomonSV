#! /bin/sh
#$ -S /bin/sh
#$ -cwd

INPUTDIR=$1
OUTPUTDIR=$2
REGION=$3

echo "samtools view ${INPUTDIR}/${REGION}.bam | perl bamPerseSV.pl - > ${OUTPUTDIR}/${REGION}.junction.temp.txt"
samtools view ${INPUTDIR}/${REGION}.bam | perl bamPerseSV.pl - > ${OUTPUTDIR}/${REGION}.junction.temp.txt

# cat ${OUTPUTDIR}/*junction.temp.txt | perl sortJunction.pl - | sort -k1,1 -k2,2n -k4,4 -k5,5n - > merge.junction.sort.txt
