#! /bin/sh
#$ -S /bin/sh
#$ -cwd

INPUTBAM=$1
OUPTUTDIR=$2
INTERVALLIST=$3

# check whether the input bam file exists or not
if [ ! -f ${INPUTBAM} ]
then
    echo "${INPUTBAM} does not exist!"
    exit
fi

# create the output directory if it does not exist
if [ ! -d ${OUTPUTDIR}/tmp ]
then
    mkdir -p ${OUPTUTDIR}/tmp
fi

REGIONCOUNT=`wc -l ${INTERVALLIST}`

# generate junction positions
job_parseJunction=junction.$(date +%s%N)
echo "qsub -t 1-${REGIONCOUNT}:1 -N ${job_parseJunction} ${LOGSTR} parseJunction.sh ${INPUTBAM} ${OUTPUTDIR} ${INTERVALLIST}"
qsub -t 1-${REGIONCOUNT}:1 -N ${job_parseJunction} ${LOGSTR} parseJunction.sh ${INPUTBAM} ${OUTPUTDIR} ${INTERVALLIST}

job_parseImproper=improper.$(date +%s%N)
echo "qsub -t 1-${REGIONCOUNT}:1 -N ${job_parseImproper} ${LOGSTR} parseImproperPair.sh ${INPUTBAM} ${OUTPUTDIR} ${INTERVALLIST}"
qsub -t 1-${REGIONCOUNT}:1 -N ${job_parseImproper} ${LOGSTR} parseImproperPair.sh ${INPUTBAM} ${OUTPUTDIR} ${INTERVALLIST}


