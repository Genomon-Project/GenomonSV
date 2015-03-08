#! /bin/sh
#$ -S /bin/sh
#$ -cwd

INPUTDIR=$1

echo "cat ${INPUTDIR}/*.juncPairInfo.txt | sort -k 5 -n > ${INPUTDIR}/merge.juncPairInfo.sort.txt"
cat ${INPUTDIR}/*.juncPairInfo.txt | sort -k 5 -n > ${INPUTDIR}/merge.juncPairInfo.sort.txt

echo "python addJuncPairInfo.py ${INPUTDIR}/merge.junction.sort.txt ${INPUTDIR}/merge.juncPairInfo.sort.txt > ${INPUTDIR}/merge.juncWithPair.txt"
python addJuncPairInfo.py ${INPUTDIR}/merge.junction.sort.txt ${INPUTDIR}/merge.juncPairInfo.sort.txt > ${INPUTDIR}/merge.juncWithPair.txt

echo "perl summarize.junctionBedpe.pl ${INPUTDIR}/merge.juncWithPair.txt > ${INPUTDIR}/merge.junction.summarized.bedpe"
perl summarize.junctionBedpe.pl ${INPUTDIR}/merge.juncWithPair.txt > ${INPUTDIR}/merge.junction.summarized.bedpe
 

