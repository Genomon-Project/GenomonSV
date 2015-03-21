#! /bin/sh
#$ -S /bin/sh
#$ -cwd

source ./config.sh

INPUTBAM=$1
OUTPUTDIR=$2
INTERVALLIST=$3
# NUM=$4

REGION=`head -n ${SGE_TASK_ID} ${INTERVALLIST} | tail -n 1`
# REGION=`head -n ${NUM} ${INTERVALLIST} | tail -n 1`

echo "python parseJunctionFromBam.py ${INPUTBAM} ${REGION} > ${OUTPUTDIR}/tmp/${REGION}.junction.temp.txt"
python parseJunctionFromBam.py ${INPUTBAM} ${REGION} > ${OUTPUTDIR}/tmp/${REGION}.junction.temp.txt

# echo "samtools view ${INPUTBAM} ${REGION} | perl bamPerseSV.pl - > ${OUTPUTDIR}/tmp/${REGION}.junction.temp.txt"
# samtools view ${INPUTBAM} ${REGION} | perl bamPerseSV.pl - > ${OUTPUTDIR}/tmp/${REGION}.junction.temp.txt

