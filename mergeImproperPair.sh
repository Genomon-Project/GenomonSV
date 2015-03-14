#! /bin/sh
#$ -S /bin/sh
#$ -cwd


OUTPUTDIR=$1

echo "cat ${OUTPUTDIR}/tmp/*improperPair.txt | sort -k1 - > ${OUTPUTDIR}/merge.improperPair.txt"
cat ${OUTPUTDIR}/tmp/*improperPair.txt | sort -k1 - > ${OUTPUTDIR}/merge.improperPair.txt

echo "perl makeBedpeFromImproperPair.pl ${OUTPUTDIR}/merge.improperPair.txt | sort -k1,1 -k2,2n -k4,4 -k5,5n - > ${OUTPUTDIR}/merge.improperPair.bedpe"
perl makeBedpeFromImproperPair.pl ${OUTPUTDIR}/merge.improperPair.txt | sort -k1,1 -k2,2n -k4,4 -k5,5n - > ${OUTPUTDIR}/merge.improperPair.bedpe

echo "perl summarize.improperPairBedpe.pl ${OUTPUTDIR}/merge.improperPair.bedpe | sort -k1,1 -k2,2n -k4,4 -k5,5n - > ${OUTPUTDIR}/merge.improperPair.summarized.bedpe"
perl summarize.improperPairBedpe.pl ${OUTPUTDIR}/merge.improperPair.bedpe | sort -k1,1 -k2,2n -k4,4 -k5,5n - > ${OUTPUTDIR}/merge.improperPair.summarized.bedpe

echo "bgzip ${OUTPUTDIR}/merge.improperPair.summarized.bedpe > ${OUTPUTDIR}/merge.improperPair.summarized.bedpe.gz"
bgzip ${OUTPUTDIR}/merge.improperPair.summarized.bedpe > ${OUTPUTDIR}/merge.improperPair.summarized.bedpe.gz

echo "tabix -p bed ${OUTPUTDIR}/merge.improperPair.summarized.bedpe.gz"
tabix -p bed ${OUTPUTDIR}/merge.improperPair.summarized.bedpe.gz



