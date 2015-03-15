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
 
# echo "tabix ${OUTPUTDIR}/merge.junctionPair.sort.bed.gz ${REGION} > ${OUTPUTDIR}/${REGION}.junctionPair.sort.bed"
# tabix ${OUTPUTDIR}/merge.junctionPair.sort.bed.gz ${REGION} > ${OUTPUTDIR}/${REGION}.junctionPair.sort.bed 

# echo "samtools view ${SEQDIR}/${REGION}.bam | perl getJuncPairInfoFromBam.pl - ${OUTPUTDIR}/${REGION}.junctionPair.sort.bed > ${OUTPUTDIR}/${REGION}.juncPairInfo.txt"
# samtools view ${SEQDIR}/${REGION}.bam | perl getJuncPairInfoFromBam.pl - ${OUTPUTDIR}/${REGION}.junctionPair.sort.bed > ${OUTPUTDIR}/${REGION}.juncPairInfo.txt

