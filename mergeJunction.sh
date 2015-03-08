#! /bin/sh
#$ -S /bin/sh
#$ -cwd

OUTPUTDIR=$1

cat ${OUTPUTDIR}/*junction.temp.txt | perl sortJunction.pl - | sort -k1,1 -k2,2n -k4,4 -k5,5n - > ${OUTPUTDIR}/merge.junction.sort.txt


echo "perl getJuncPair.pl ${OUTPUTDIR}/merge.junction.sort.txt | sort -k1,1 -k3,3n > ${OUTPUTDIR}/merge.junctionPair.sort.bed"
perl getJuncPair.pl ${OUTPUTDIR}/merge.junction.sort.txt | sort -k1,1 -k3,3n > ${OUTPUTDIR}/merge.junctionPair.sort.bed

echo "bgzip ${OUTPUTDIR}/merge.junctionPair.sort.bed > ${OUTPUTDIR}/merge.junctionPair.sort.bed.gz"
bgzip ${OUTPUTDIR}/merge.junctionPair.sort.bed > ${OUTPUTDIR}/merge.junctionPair.sort.bed.gz


tabix -p bed ${OUTPUTDIR}/merge.junctionPair.sort.bed.gz

