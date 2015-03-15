#! /bin/sh
#$ -S /bin/sh
#$ -cwd

export PATH=/home/yshira/bin/tabix-0.2.6:$PATH

SEQDIR=$1
OUTPUTDIR=$2
REGION=$3

echo "python getPairInfoFromBam.py ${INPUTBAM} ${OUTPUTDIR}/merge.junctionPair.sort.bed.gz ${REGION} > ${OUTPUTDIR}/${REGION}.juncPairInfo.txt"
python getPairInfoFromBam.py ${INPUTBAM} ${OUTPUTDIR}/merge.junctionPair.sort.bed.gz ${REGION} > ${OUTPUTDIR}/${REGION}.juncPairInfo.txt
 
# echo "tabix ${OUTPUTDIR}/merge.junctionPair.sort.bed.gz ${REGION} > ${OUTPUTDIR}/${REGION}.junctionPair.sort.bed"
# tabix ${OUTPUTDIR}/merge.junctionPair.sort.bed.gz ${REGION} > ${OUTPUTDIR}/${REGION}.junctionPair.sort.bed 

# echo "samtools view ${SEQDIR}/${REGION}.bam | perl getJuncPairInfoFromBam.pl - ${OUTPUTDIR}/${REGION}.junctionPair.sort.bed > ${OUTPUTDIR}/${REGION}.juncPairInfo.txt"
# samtools view ${SEQDIR}/${REGION}.bam | perl getJuncPairInfoFromBam.pl - ${OUTPUTDIR}/${REGION}.junctionPair.sort.bed > ${OUTPUTDIR}/${REGION}.juncPairInfo.txt

