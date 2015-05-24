#! /bin/sh
#$ -S /bin/sh
#$ -cwd

##########
# The script performs
# 1. extract the read pairs around the break point
# 2. make the reference sequence and the sequence with the structural variation candidate
# 3. alignment the short read pairs by blat
# 4. summarized the aligned read pairs
##########

SUFFIX=`sh getSuffix.sh ${SGE_TASK_ID}`

INPUT=$1.${SUFFIX}
TUMORBAM=$2
NORMALBAM=$3
OUTPUT=$4.${SUFFIX}

source ./config.sh

echo -n > ${OUTPUT}
while read LINE; do 

    chr1=`echo ${LINE} | cut -d ' ' -f 1`
    pos1=`echo ${LINE} | cut -d ' ' -f 3`
    dir1=`echo ${LINE} | cut -d ' ' -f 9`
    chr2=`echo ${LINE} | cut -d ' ' -f 4`
    pos2=`echo ${LINE} | cut -d ' ' -f 6`
    dir2=`echo ${LINE} | cut -d ' ' -f 10`
    juncSeq=`echo ${LINE} | cut -d ' ' -f 8`

    ITDFlag=0
    if [ $chr1 = $chr2 -a `expr $pos2 - $pos1` -lt 500 -a $dir1 = "-" -a $dir2 = "+" ];
    then
        ITDFlag=1
    fi

    echo "python extractSVReadPairs.py ${TUMORBAM} ${chr1} ${pos1} ${dir1} ${chr2} ${pos2} ${dir2} 1000 5 > ${OUTPUT}.temp.tumor.fa"
    python extractSVReadPairs.py ${TUMORBAM} ${chr1} ${pos1} ${dir1} ${chr2} ${pos2} ${dir2} 1000 5 > ${OUTPUT}.temp.tumor.fa
    if [ $? -eq  27 ]; then
        echo -e "${chr1}\t${pos1}\t${dir1}\t${chr2}\t${pos2}\t${dir2}\t${juncSeq}\t*\t*\t*\t*\t*" >> ${OUTPUT}
        continue
    fi
    check_error $?

    echo "python extractSVReadPairs.py ${NORMALBAM} ${chr1} ${pos1} ${dir1} ${chr2} ${pos2} ${dir2} 1000 5 > ${OUTPUT}.temp.normal.fa"
    python extractSVReadPairs.py ${NORMALBAM} ${chr1} ${pos1} ${dir1} ${chr2} ${pos2} ${dir2} 1000 5 > ${OUTPUT}.temp.normal.fa
    if [ $? -eq  27 ]; then
        echo -e "${chr1}\t${pos1}\t${dir1}\t${chr2}\t${pos2}\t${dir2}\t${juncSeq}\t*\t*\t*\t*\t*" >> ${OUTPUT}
        continue
    fi
    check_error $?


    echo "python getRefAltForSV.py ${REFERENCE} ${chr1} ${pos1} ${dir1} ${chr2} ${pos2} ${dir2} ${juncSeq} > ${OUTPUT}.temp.refalt.fa"
    python getRefAltForSV.py ${REFERENCE} ${chr1} ${pos1} ${dir1} ${chr2} ${pos2} ${dir2} ${juncSeq} > ${OUTPUT}.temp.refalt.fa
    check_error $?

    # wc ${OUTPUT}.temp.tumor.fa
    echo "blat -stepSize=5 -repMatch=2253 ${OUTPUT}.temp.refalt.fa ${OUTPUT}.temp.tumor.fa ${OUTPUT}.temp.tumor.psl"
    blat -stepSize=5 -repMatch=2253 ${OUTPUT}.temp.refalt.fa ${OUTPUT}.temp.tumor.fa ${OUTPUT}.temp.tumor.psl
    check_error $?

    # wc ${OUTPUT}.temp.normal.fa
    echo "blat -stepSize=5 -repMatch=2253 ${OUTPUT}.temp.refalt.fa ${OUTPUT}.temp.normal.fa ${OUTPUT}.temp.normal.psl"
    blat -stepSize=5 -repMatch=2253 ${OUTPUT}.temp.refalt.fa ${OUTPUT}.temp.normal.fa ${OUTPUT}.temp.normal.psl
    check_error $?

    echo "python procPslFisher.py ${OUTPUT}.temp.tumor.psl ${OUTPUT}.temp.normal.psl ${ITDFlag} > ${OUTPUT}.temp.fisher.txt"
    python procPslFisher.py ${OUTPUT}.temp.tumor.psl ${OUTPUT}.temp.normal.psl ${ITDFlag} > ${OUTPUT}.temp.fisher.txt
    check_error $?

    echo -ne "${chr1}\t${pos1}\t${dir1}\t${chr2}\t${pos2}\t${dir2}\t${juncSeq}\t" >> ${OUTPUT}
    cat ${OUTPUT}.temp.fisher.txt >> ${OUTPUT}
    check_error $?

done < ${INPUT}


