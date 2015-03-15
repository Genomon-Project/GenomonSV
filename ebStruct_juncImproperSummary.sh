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
job_parseJunction=junction.$(date +%s%N)
echo "qsub -t 1-${REGIONCOUNT}:1 -N ${job_parseJunction} -e log/ -o log/ parseJunction.sh ${INPUTBAM} ${OUTPUTDIR} ${INTERVALLIST}"
qsub -t 1-${REGIONCOUNT}:1 -N ${job_parseJunction} -e log/ -o log/ parseJunction.sh ${INPUTBAM} ${OUTPUTDIR} ${INTERVALLIST}

job_parseImproper=improper.$(date +%s%N)
echo "qsub -t 1-${REGIONCOUNT}:1 -N ${job_parseImproper} -e log/ -o log/ parseImproperPair.sh ${INPUTBAM} ${OUTPUTDIR} ${INTERVALLIST}"
qsub -t 1-${REGIONCOUNT}:1 -N ${job_parseImproper} -e log/ -o log/ parseImproperPair.sh ${INPUTBAM} ${OUTPUTDIR} ${INTERVALLIST}


