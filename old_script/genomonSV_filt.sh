#! /bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -e log/ -o log/


TUMORDIR=$1
TUMORBAM=$2
NORMALBAM=$3
CONTROL=$4
MATCHEDNORMAL=$5

# check whether the input files exist or not
if [ ! -f ${TUMORDIR}/merge.junction.summarized.bedpe.gz ]
then
    echo "${TUMORDIR}/merge.junction.summarized.bedpe.gz does not exist!"
    exit
fi

if [ ! -f ${TUMORBAM} ]
then
    echo "${TUMORBAM} does not exist!"
    exit
fi

if [ ! -f ${NORMALBAM} ]
then 
    echo "${NORMALBAM} does not exist!"
    exit
fi  

if [ ! -f ${CONTROL} ]
then
    echo "${CONTROL} does not exist!"
    exit
fi


# generate somatic SV candidate 
job_filterCandidate=filterCandidate.$(date +%s%N)
echo "qsub -e log/ -o log/ -N ${job_filterCandidate} filterCandidate.sh ${TUMORDIR} ${CONTROL} ${MATCHEDNORMAL}"
qsub -e log/ -o log/ -N ${job_filterCandidate} filterCandidate.sh ${TUMORDIR} ${CONTROL} ${MATCHEDNORMAL}

job_master_validate=master_validate.$(date +%s%N)
echo "qsub -e log/ -o log/ -N ${job_master_validate} -hold_jid ${job_filterCandidate} master_validateByRealignment.sh ${TUMORDIR} ${TUMORBAM} ${NORMALBAM}"
qsub -e log/ -o log/ -N ${job_master_validate} -hold_jid ${job_filterCandidate} master_validateByRealignment.sh ${TUMORDIR} ${TUMORBAM} ${NORMALBAM}


# job_anno=annotAndFilter.$(date +%s%N)
# echo "qsub -e log/ -o log/ -N ${job_anno} catCandSV.sh ${TUMORDIR}"
# qsub -e log/ -o log/ -N ${job_anno} catCandSV.sh ${TUMORDIR}

