#! /bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -e log/ -o log/

source ./config.sh

TUMORDIR=$1
NORMALDIR=$2
TUMORBAM=$3
NORMALBAM=$4
CONTROL=$5
MATCHEDNORMAL=$6

echo "python filterLengthNum.py ${TUMORDIR}/merge.junction.summarized.bedpe.gz 2 1000 > ${TUMORDIR}/merge.junction.summarized.filt1.bedpe"
python filterLengthNum.py ${TUMORDIR}/merge.junction.summarized.bedpe.gz 2 1000 > ${TUMORDIR}/merge.junction.summarized.filt1.bedpe
check_error $?

echo "python filterNonMatchControl.py ${TUMORDIR}/merge.junction.summarized.filt1.bedpe ${CONTROL} ${MATCHEDNORMAL} 1 > ${TUMORDIR}/merge.junction.summarized.filt2.bedpe"
python filterNonMatchControl.py ${TUMORDIR}/merge.junction.summarized.filt1.bedpe ${CONTROL} ${MATCHEDNORMAL} 1 > ${TUMORDIR}/merge.junction.summarized.filt2.bedpe 
check_error $?

echo "python addImproperInfo.py ${TUMORDIR}/merge.junction.summarized.filt2.bedpe ${TUMORDIR}/merge.improperPair.summarized.bedpe.gz > ${TUMORDIR}/merge.junction.summarized.filt3.bedpe"
python addImproperInfo.py ${TUMORDIR}/merge.junction.summarized.filt2.bedpe ${TUMORDIR}/merge.improperPair.summarized.bedpe.gz > ${TUMORDIR}/merge.junction.summarized.filt3.bedpe 
check_error $?

echo "python filterMergedJunc.py ${TUMORDIR}/merge.junction.summarized.filt3.bedpe 3 40 100 > ${TUMORDIR}/merge.junction.summarized.filt4.bedpe"
python filterMergedJunc.py ${TUMORDIR}/merge.junction.summarized.filt3.bedpe 3 40 100 > ${TUMORDIR}/merge.junction.summarized.filt4.bedpe 
check_error $?

echo "python filterFisher.py ${TUMORDIR}/merge.junction.summarized.filt4.bedpe ${NORMALDIR}/merge.junction.summarized.bedpe.gz ${NORMALDIR}/merge.improperPair.summarized.bedpe.gz ${TUMORBAM} ${NORMALBAM} > ${TUMORDIR}/merge.junction.summarized.filt5.bedpe"
python filterFisher.py ${TUMORDIR}/merge.junction.summarized.filt4.bedpe ${NORMALDIR}/merge.junction.summarized.bedpe.gz ${NORMALDIR}/merge.improperPair.summarized.bedpe.gz ${TUMORBAM} ${NORMALBAM} > ${TUMORDIR}/merge.junction.summarized.filt5.bedpe
check_error $?

