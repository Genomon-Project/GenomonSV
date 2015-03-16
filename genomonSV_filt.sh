#! /bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -e log/ -o log/

source ./config.sh

OUTPUTDIR=$1
CONTROL=$2
MATCHEDNORMAL=$3

echo "python filterLengthNum.py ${OUTPUTDIR}/merge.junction.summarized.bedpe 2 1000 > ${OUTPUTDIR}/merge.junction.summarized.filt1.bedpe"
python filterLengthNum.py ${OUTPUTDIR}/merge.junction.summarized.bedpe 2 1000 > ${OUTPUTDIR}/merge.junction.summarized.filt1.bedpe

echo "python filterNonMatchControl.py ${OUTPUTDIR}/merge.junction.summarized.filt1.bedpe ${CONTROL} ${MATCHEDNORMAL} 1 > ${OUTPUTDIR}/merge.junction.summarized.filt2.bedpe"
python filterNonMatchControl.py ${OUTPUTDIR}/merge.junction.summarized.filt1.bedpe ${CONTROL} ${MATCHEDNORMAL} 1 > ${OUTPUTDIR}/merge.junction.summarized.filt2.bedpe 

echo "python addImproperInfo.py ${OUTPUTDIR}/merge.junction.summarized.filt2.bedpe ${OUTPUTDIR}/merge.improperPair.summarized.bedpe.gz > ${OUTPUTDIR}/merge.junction.summarized.filt3.bedpe"
python addImproperInfo.py ${OUTPUTDIR}/merge.junction.summarized.filt2.bedpe ${OUTPUTDIR}/merge.improperPair.summarized.bedpe.gz > ${OUTPUTDIR}/merge.junction.summarized.filt3.bedpe 

echo "python filterMergedJunc.py ${OUTPUTDIR}/merge.junction.summarized.filt3.bedpe 3 40 100 > ${OUTPUTDIR}/merge.junction.summarized.filt4.bedpe"
python filterMergedJunc.py ${OUTPUTDIR}/merge.junction.summarized.filt3.bedpe 3 40 100 > ${OUTPUTDIR}/merge.junction.summarized.filt4.bedpe 

