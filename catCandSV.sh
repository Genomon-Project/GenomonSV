#! /bin/sh
#$ -S /bin/sh
#$ -cwd

TUMORDIR=$1

source ./config.sh

echo "cat ${TUMORDIR}/fisherTmp/svCand.fisher.txt.??? > ${TUMORDIR}/fisherTmp/svCand.fisher.txt"
cat ${TUMORDIR}/fisherTmp/svCand.fisher.txt.??? > ${TUMORDIR}/fisherTmp/svCand.fisher.txt
check_error $?

echo "python filterAndAnno.py ${TUMORDIR}/fisherTmp/svCand.fisher.txt db/refGene.bed.gz db/refExon.bed.gz > ${TUMORDIR}/genomonSV.result.txt"
python filterAndAnno.py ${TUMORDIR}/fisherTmp/svCand.fisher.txt db/refGene.bed.gz db/refExon.bed.gz > ${TUMORDIR}/genomonSV.result.txt
check_error $?


