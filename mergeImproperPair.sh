#! /bin/sh
#$ -S /bin/sh
#$ -cwd

source ./config.sh

OUTPUTDIR=$1

echo "cat ${OUTPUTDIR}/tmp/*improperPair.temp.txt | sort -k1 - > ${OUTPUTDIR}/merge.improperPair.txt"
cat ${OUTPUTDIR}/tmp/*improperPair.temp.txt | sort -k1 - > ${OUTPUTDIR}/merge.improperPair.txt

echo "python makeBedpeFromImproperPair.py ${OUTPUTDIR}/merge.improperPair.txt | sort -k1,1 -k2,2n -k4,4 -k5,5n - > ${OUTPUTDIR}/merge.improperPair.bedpe"
python makeBedpeFromImproperPair.py ${OUTPUTDIR}/merge.improperPair.txt | sort -k1,1 -k2,2n -k4,4 -k5,5n - > ${OUTPUTDIR}/merge.improperPair.bedpe

echo "perl summarize.improperPairBedpe.pl ${OUTPUTDIR}/merge.improperPair.bedpe | sort -k1,1 -k2,2n -k4,4 -k5,5n - > ${OUTPUTDIR}/merge.improperPair.summarized.bedpe"
perl summarize.improperPairBedpe.pl ${OUTPUTDIR}/merge.improperPair.bedpe | sort -k1,1 -k2,2n -k4,4 -k5,5n - > ${OUTPUTDIR}/merge.improperPair.summarized.bedpe

echo "bgzip -f ${OUTPUTDIR}/merge.improperPair.summarized.bedpe > ${OUTPUTDIR}/merge.improperPair.summarized.bedpe.gz"
bgzip -f ${OUTPUTDIR}/merge.improperPair.summarized.bedpe > ${OUTPUTDIR}/merge.improperPair.summarized.bedpe.gz

echo "tabix -p bed ${OUTPUTDIR}/merge.improperPair.summarized.bedpe.gz"
tabix -p bed ${OUTPUTDIR}/merge.improperPair.summarized.bedpe.gz



