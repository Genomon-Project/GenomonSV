#! /bin/sh
#$ -S /bin/sh
#$ -cwd

INPUTBAM=$1
OUTPUTDIR=$2
INTERVALLIST=$3

REGION=`head -n ${SGE_TASK_ID} ${INTERVALLIST} | tail -n 1`

echo "samtools view ${INPUTBAM} ${REGION} | perl bamPerseSV.pl - > ${OUTPUTDIR}/temp/${REGION}.junction.temp.txt"
samtools view ${INPUTBAM} ${REGION} | perl bamPerseSV.pl - > ${OUTPUTDIR}/temp/${REGION}.junction.temp.txt

