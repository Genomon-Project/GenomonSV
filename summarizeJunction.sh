#! /bin/sh
#$ -S /bin/sh
#$ -cwd

source ./config.sh

INPUTDIR=$1

echo "cat ${INPUTDIR}/tmp/*.juncPairInfo.txt | sort -k 5 -n > ${INPUTDIR}/merge.juncPairInfo.sort.txt"
cat ${INPUTDIR}/tmp/*.juncPairInfo.txt | sort -k 5 -n > ${INPUTDIR}/merge.juncPairInfo.sort.txt
check_error $?

echo "python addJuncPairInfo.py ${INPUTDIR}/merge.junction.sort.txt ${INPUTDIR}/merge.juncPairInfo.sort.txt > ${INPUTDIR}/merge.juncWithPair.txt"
python addJuncPairInfo.py ${INPUTDIR}/merge.junction.sort.txt ${INPUTDIR}/merge.juncPairInfo.sort.txt > ${INPUTDIR}/merge.juncWithPair.txt
check_error $?

echo "python summarizeJunctionBedpe.py ${INPUTDIR}/merge.juncWithPair.txt | sort -k1,1 -k2,2n -k4,4 -k5,5n - > ${INPUTDIR}/merge.junction.summarized.bedpe"
python summarizeJunctionBedpe.py ${INPUTDIR}/merge.juncWithPair.txt | sort -k1,1 -k2,2n -k4,4 -k5,5n - > ${INPUTDIR}/merge.junction.summarized.bedpe
check_error $?

echo "bgzip -f ${INPUTDIR}/merge.junction.summarized.bedpe > ${INPUTDIR}/merge.junction.summarized.bedpe.gz" 
bgzip -f ${INPUTDIR}/merge.junction.summarized.bedpe > ${INPUTDIR}/merge.junction.summarized.bedpe.gz
check_error $?

echo "tabix -p bed ${INPUTDIR}/merge.junction.summarized.bedpe.gz"
tabix -p bed ${INPUTDIR}/merge.junction.summarized.bedpe.gz
check_error $?


