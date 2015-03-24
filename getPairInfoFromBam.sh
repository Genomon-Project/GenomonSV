#! /bin/sh
#$ -S /bin/sh
#$ -cwd

source ./config.sh

INPUTBAM=$1
OUTPUTDIR=$2
INTERVALLIST=$3
REGION=`head -n ${SGE_TASK_ID} ${INTERVALLIST} | tail -n 1`

echo "python getPairInfoFromBam.py ${INPUTBAM} ${OUTPUTDIR}/merge.junctionPair.sort.bed.gz ${REGION} > ${OUTPUTDIR}/tmp/${REGION}.juncPairInfo.txt"
python getPairInfoFromBam.py ${INPUTBAM} ${OUTPUTDIR}/merge.junctionPair.sort.bed.gz ${REGION} > ${OUTPUTDIR}/tmp/${REGION}.juncPairInfo.txt
check_error $?
 

