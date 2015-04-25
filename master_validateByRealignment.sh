#! /bin/sh
#$ -S /bin/sh
#$ -cwd

TUMORDIR=$1
TUMORBAM=$2
NORMALBAM=$3

INPUTCOUNT=`ls ${TUMORDIR}/fisherTmp/svCand.txt.* | wc -l`

job_validateByRealignment=validateByRealignmentdate.$(date +%s%N)
echo "qsub -sync y -t 1-${INPUTCOUNT}:1 -N ${job_validateByRealignment} -e log/ -o log/ validateByRealignment.sh ${TUMORDIR}/fisherTmp/svCand.txt ${TUMORBAM} ${NORMALBAM} ${TUMORDIR}/fisherTmp/svCand.fisher.txt"
qsub -sync y -t 1-${INPUTCOUNT}:1 -N ${job_validateByRealignment} -e log/ -o log/ validateByRealignment.sh ${TUMORDIR}/fisherTmp/svCand.txt ${TUMORBAM} ${NORMALBAM} ${TUMORDIR}/fisherTmp/svCand.fisher.txt

