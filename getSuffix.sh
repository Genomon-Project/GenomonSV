#! /bin/sh
#$ -S /bin/sh
#$ -cwd

NUM=$1
NUM=`expr ${NUM} - 1`

array=( a b c d e f g h i j k l m n o p q r s t u v w x y z )

DIGIT1=`expr ${NUM} / 676`
NUM1=`expr ${NUM} % 676`

DIGIT2=`expr ${NUM1} / 26`
DIGIT3=`expr ${NUM1} % 26`

echo ${array[$DIGIT1]}${array[$DIGIT2]}${array[$DIGIT3]} 

