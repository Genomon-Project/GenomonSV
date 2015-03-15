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

# echo "samtools view ${INPUTBAM} ${REGION} | perl bamPerseSV2.pl - > ${OUTPUTDIR}/${REGION}.improperPair.txt"
# samtools view ${INPUTBAM} ${REGION} | perl bamPerseSV2.pl - > ${OUTPUTDIR}/${REGION}.improperPair.txt

echo "python bamParseSV2.py ${INPUTBAM} ${REGION} > ${OUTPUTDIR}/tmp/${REGION}.improperPair.temp.txt"
python bamParseSV2.py ${INPUTBAM} ${REGION} > ${OUTPUTDIR}/tmp/${REGION}.improperPair.temp.txt
 
