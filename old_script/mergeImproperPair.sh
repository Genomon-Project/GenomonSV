#! /bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -q mjobs.q

source ./config.sh

OUTPUTDIR=$1

echo "cat ${OUTPUTDIR}/tmp/*improperPair.temp.txt | sort -k1 - > ${OUTPUTDIR}/merge.improperPair.txt"
cat ${OUTPUTDIR}/tmp/*improperPair.temp.txt | sort -k1 - > ${OUTPUTDIR}/merge.improperPair.txt
check_error $?

echo "python makeBedpeFromImproperPair.py ${OUTPUTDIR}/merge.improperPair.txt | sort -k1,1 -k2,2n -k4,4 -k5,5n - > ${OUTPUTDIR}/merge.improperPair.bedpe"
python makeBedpeFromImproperPair.py ${OUTPUTDIR}/merge.improperPair.txt | sort -k1,1 -k2,2n -k4,4 -k5,5n - > ${OUTPUTDIR}/merge.improperPair.bedpe
check_error $?

echo "python summarizeImproperPairBedpe.py ${OUTPUTDIR}/merge.improperPair.bedpe | sort -k1,1 -k2,2n -k4,4 -k5,5n - > ${OUTPUTDIR}/merge.improperPair.summarized.bedpe"
python summarizeImproperPairBedpe.py ${OUTPUTDIR}/merge.improperPair.bedpe | sort -k1,1 -k2,2n -k4,4 -k5,5n - > ${OUTPUTDIR}/merge.improperPair.summarized.bedpe
check_error $?

echo "bgzip -f ${OUTPUTDIR}/merge.improperPair.summarized.bedpe > ${OUTPUTDIR}/merge.improperPair.summarized.bedpe.gz"
bgzip -f ${OUTPUTDIR}/merge.improperPair.summarized.bedpe > ${OUTPUTDIR}/merge.improperPair.summarized.bedpe.gz
check_error $?

echo "tabix -p bed ${OUTPUTDIR}/merge.improperPair.summarized.bedpe.gz"
tabix -p bed ${OUTPUTDIR}/merge.improperPair.summarized.bedpe.gz
check_error $?



