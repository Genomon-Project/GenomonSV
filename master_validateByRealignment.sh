#! /bin/sh
#$ -S /bin/sh
#$ -cwd

TUMORDIR=$1
TUMORBAM=$2
NORMALBAM=$3

INPUTCOUNT=`ls ${TUMORDIR}/fisherTmp/svCand.txt.* | wc -l`

if [ ${INPUTCOUNT} -gt 0 ];
then
    job_validateByRealignment=validateByRealignmentdate.$(date +%s%N)
    echo "qsub -t 1-${INPUTCOUNT}:1 -N ${job_validateByRealignment} -e log/ -o log/ validateByRealignment.sh ${TUMORDIR}/fisherTmp/svCand.txt ${TUMORBAM} ${NORMALBAM} ${TUMORDIR}/fisherTmp/svCand.fisher.txt"
    qsub -t 1-${INPUTCOUNT}:1 -N ${job_validateByRealignment} -e log/ -o log/ validateByRealignment.sh ${TUMORDIR}/fisherTmp/svCand.txt ${TUMORBAM} ${NORMALBAM} ${TUMORDIR}/fisherTmp/svCand.fisher.txt

    job_anno=annotAndFilter.$(date +%s%N)
    echo "qsub -e log/ -o log/ -N ${job_anno} -hold_jid ${job_validateByRealignment} catCandSV.sh ${TUMORDIR}"
    qsub -e log/ -o log/ -N ${job_anno} -hold_jid ${job_validateByRealignment} catCandSV.sh ${TUMORDIR}

fi

