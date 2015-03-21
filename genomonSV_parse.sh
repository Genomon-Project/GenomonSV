#! /bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -e log/ -o log/

INPUTBAM=$1
OUTPUTDIR=$2
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
    mkdir -p ${OUTPUTDIR}/tmp
fi

REGIONCOUNT=`wc -l ${INTERVALLIST} | cut -d ' ' -f 1`

# generate junction positions
job_parseJunction=parseJunction.$(date +%s%N)
echo "qsub -t 1-${REGIONCOUNT}:1 -N ${job_parseJunction} -e log/ -o log/ parseJunction.sh ${INPUTBAM} ${OUTPUTDIR} ${INTERVALLIST}"
qsub -t 1-${REGIONCOUNT}:1 -N ${job_parseJunction} -e log/ -o log/ parseJunction.sh ${INPUTBAM} ${OUTPUTDIR} ${INTERVALLIST}

job_parseImproper=parseImproper.$(date +%s%N)
echo "qsub -t 1-${REGIONCOUNT}:1 -N ${job_parseImproper} -e log/ -o log/ parseImproperPair.sh ${INPUTBAM} ${OUTPUTDIR} ${INTERVALLIST}"
qsub -t 1-${REGIONCOUNT}:1 -N ${job_parseImproper} -e log/ -o log/ parseImproperPair.sh ${INPUTBAM} ${OUTPUTDIR} ${INTERVALLIST}

job_mergeJunction=mergeJunction.$(date +%s%N)
echo "qsub -N ${job_mergeJunction} -hold_jid ${job_parseJunction} -e log/ -o log/ mergeJunction.sh ${OUTPUTDIR}"
qsub -N ${job_mergeJunction} -hold_jid ${job_parseJunction} -e log/ -o log/ mergeJunction.sh ${OUTPUTDIR}

job_mergeImproper=mergeImproper.$(date +%s%N)
echo "qsub -N ${job_mergeImproper} -hold_jid ${job_parseImproper} -e log/ -o log/ mergeImproperPair.sh ${OUTPUTDIR}"
qsub -N ${job_mergeImproper} -hold_jid ${job_parseImproper} -e log/ -o log/ mergeImproperPair.sh ${OUTPUTDIR}

job_getPairInfoJunc=getPairInfoJunc.$(date +%s%N)
echo "qsub -t 1-${REGIONCOUNT}:1 -N ${job_getPairInfoJunc} -hold_jid ${job_mergeJunction} -e log/ -o log/ getPairInfoFromBam.sh ${INPUTBAM} ${OUTPUTDIR} ${INTERVALLIST}"
# qsub -t 1-${REGIONCOUNT}:1 -N ${job_getPairInfoJunc} -hold_jid ${job_mergeJunction} -e log/ -o log/ getPairInfoFromBam.sh ${INPUTBAM} ${OUTPUTDIR} ${INTERVALLIST}

job_summarizeJunction=summarizeJunction.$(date +%s%N)
echo "qsub -N ${job_summarizeJunction} -hold_jid ${job_mergeJunction},${job_mergeImproper} -e log/ -o log/ summarizeJunction.sh ${OUTPUTDIR}"
# qsub -N ${job_summarizeJunction} -hold_jid ${job_mergeJunction},${job_mergeImproper} -e log/ -o log/ summarizeJunction.sh ${OUTPUTDIR}


