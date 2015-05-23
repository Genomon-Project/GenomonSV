#! /bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -e log/ -o log/
#$ -q mjobs.q

source ./config.sh

TUMORDIR=$1
CONTROL=$2
MATCHEDNORMAL=$3

echo "python filterLengthNum.py ${TUMORDIR}/merge.junction.summarized.bedpe.gz 5 10 > ${TUMORDIR}/merge.junction.summarized.filt1.bedpe"
python filterLengthNum.py ${TUMORDIR}/merge.junction.summarized.bedpe.gz 5 10 > ${TUMORDIR}/merge.junction.summarized.filt1.bedpe
check_error $?

echo "python filterNonMatchControl.py ${TUMORDIR}/merge.junction.summarized.filt1.bedpe ${CONTROL} ${MATCHEDNORMAL} 3 > ${TUMORDIR}/merge.junction.summarized.filt2.bedpe"
python filterNonMatchControl.py ${TUMORDIR}/merge.junction.summarized.filt1.bedpe ${CONTROL} ${MATCHEDNORMAL} 3 > ${TUMORDIR}/merge.junction.summarized.filt2.bedpe 
check_error $?

echo "python addImproperInfo.py ${TUMORDIR}/merge.junction.summarized.filt2.bedpe ${TUMORDIR}/merge.improperPair.summarized.bedpe.gz > ${TUMORDIR}/merge.junction.summarized.filt3.bedpe"
python addImproperInfo.py ${TUMORDIR}/merge.junction.summarized.filt2.bedpe ${TUMORDIR}/merge.improperPair.summarized.bedpe.gz > ${TUMORDIR}/merge.junction.summarized.filt3.bedpe 
check_error $?

echo "python filterMergedJunc.py ${TUMORDIR}/merge.junction.summarized.filt3.bedpe 5 40 100 > ${TUMORDIR}/merge.junction.summarized.filt4.bedpe"
python filterMergedJunc.py ${TUMORDIR}/merge.junction.summarized.filt3.bedpe 5 40 100 > ${TUMORDIR}/merge.junction.summarized.filt4.bedpe 
check_error $?

echo "python removeClose.py ${TUMORDIR}/merge.junction.summarized.filt4.bedpe > ${TUMORDIR}/merge.junction.summarized.filt5.bedpe"
python removeClose.py ${TUMORDIR}/merge.junction.summarized.filt4.bedpe > ${TUMORDIR}/merge.junction.summarized.filt5.bedpe
check_error $?

echo "rm -rf ${TUMORDIR}/fisherTmp"
rm -rf ${TUMORDIR}/fisherTmp

echo "mkdir ${TUMORDIR}/fisherTmp"
mkdir ${TUMORDIR}/fisherTmp

echo "split -l 300 -a 3 ${SUFFIX} ${TUMORDIR}/merge.junction.summarized.filt5.bedpe ${TUMORDIR}/fisherTmp/svCand.txt."
split -l 300 -a 3 ${SUFFIX} ${TUMORDIR}/merge.junction.summarized.filt5.bedpe ${TUMORDIR}/fisherTmp/svCand.txt.


# echo "python filterFisher.py ${TUMORDIR}/merge.junction.summarized.filt4.bedpe ${NORMALDIR}/merge.junction.summarized.bedpe.gz ${NORMALDIR}/merge.improperPair.summarized.bedpe.gz ${TUMORBAM} ${NORMALBAM} > ${TUMORDIR}/merge.junction.summarized.filt5.bedpe"
# python filterFisher.py ${TUMORDIR}/merge.junction.summarized.filt4.bedpe ${NORMALDIR}/merge.junction.summarized.bedpe.gz ${NORMALDIR}/merge.improperPair.summarized.bedpe.gz ${TUMORBAM} ${NORMALBAM} > ${TUMORDIR}/merge.junction.summarized.filt5.bedpe
# check_error $?

