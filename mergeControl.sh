#! /bin/sh
#$ -S /bin/sh
#$ -cwd

:<<_COMMENT_OUT_
    script for creating merged nonmatched normal control data
    the input ${CONTROLLIST} is tab-delimited file where the 1st colum is the label of samples
    and the 2nd column is the path for the individual junction list file
_COMMENT_OUT_

CONTROLLIST=$1
OUTPUT=$2

# :<<_COMMENT_OUT_

IFS=$'\t'
echo -n > ${OUTPUT}.temp
while read line
do
    tsvList=(`echo "$line"`)
    sample=${tsvList[0]}
    path=${tsvList[1]}

    echo $sample

    echo "python simplifyJunc.py $path $sample >> ${OUTPUT}.temp"
    python simplifyJunc.py $path $sample >> ${OUTPUT}.temp
 
done < ${CONTROLLIST}

echo "sort -k1,1 -k2,2n -k4,4 -k5,5n ${OUTPUT}.temp > ${OUTPUT}.temp.sort"
sort -k1,1 -k2,2n -k4,4 -k5,5n ${OUTPUT}.temp > ${OUTPUT}.temp.sort

# _COMMENT_OUT_

echo "python organizeControl.py ${OUTPUT}.temp.sort > ${OUTPUT}"
python organizeControl.py ${OUTPUT}.temp.sort > ${OUTPUT}

echo "bgzip -f ${OUTPUT} > ${OUTPUT}.gz"
bgzip -f ${OUTPUT} > ${OUTPUT}.gz

echo "tabix -p bed -f ${OUTPUT}.gz"
tabix -p bed -f ${OUTPUT}.gz


rm -rf ${OUTPUT}.temp
rm -rf ${OUTPUT}.temp.sort

