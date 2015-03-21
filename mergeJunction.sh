#! /bin/sh
#$ -S /bin/sh
#$ -cwd

source ./config.sh

OUTPUTDIR=$1

# echo "cat ${OUTPUTDIR}/tmp/*junction.temp.txt > ${OUTPUTDIR}/tmp/merge.junction.temp.txt"
# cat ${OUTPUTDIR}/tmp/*junction.temp.txt > ${OUTPUTDIR}/tmp/merge.junction.temp.txt 

# echo "python reorderJunction.py ${OUTPUTDIR}/tmp/merge.junction.temp.txt | sort -k1,1 -k2,2n -k4,4 -k5,5n - > ${OUTPUTDIR}/merge.junction.sort.txt"
# python reorderJunction.py ${OUTPUTDIR}/tmp/merge.junction.temp.txt | sort -k1,1 -k2,2n -k4,4 -k5,5n - > ${OUTPUTDIR}/merge.junction.sort.txt

echo "cat ${OUTPUTDIR}/tmp/*junction.temp.txt | sort -k1,1 -k2,2n -k4,4 -k5,5n - > ${OUTPUTDIR}/merge.junction.sort.txt"
cat ${OUTPUTDIR}/tmp/*junction.temp.txt | sort -k1,1 -k2,2n -k4,4 -k5,5n - > ${OUTPUTDIR}/merge.junction.sort.txt

echo "python getPairStartPos.py ${OUTPUTDIR}/merge.junction.sort.txt | sort -k1,1 -k3,3n > ${OUTPUTDIR}/merge.junctionPair.sort.bed"
python getPairStartPos.py ${OUTPUTDIR}/merge.junction.sort.txt | sort -k1,1 -k3,3n > ${OUTPUTDIR}/merge.junctionPair.sort.bed

echo "bgzip -f ${OUTPUTDIR}/merge.junctionPair.sort.bed > ${OUTPUTDIR}/merge.junctionPair.sort.bed.gz"
bgzip -f ${OUTPUTDIR}/merge.junctionPair.sort.bed > ${OUTPUTDIR}/merge.junctionPair.sort.bed.gz

echo "tabix -p bed ${OUTPUTDIR}/merge.junctionPair.sort.bed.gz"
tabix -p bed ${OUTPUTDIR}/merge.junctionPair.sort.bed.gz

# rm -rf ${OUTPUTDIR}/tmp/merge.junction.temp.txt

