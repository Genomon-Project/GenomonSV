#! /bin/sh


wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz

echo "python listRefGene.py refGene.txt.gz | sort -k1,1 -k2,2n -k3,3n - > refGene.bed"
python listRefGene.py refGene.txt.gz | sort -k1,1 -k2,2n -k3,3n - > refGene.bed

echo "python listRefExon.py refGene.txt.gz | sort -k1,1 -k2,2n -k3,3n - > refExon.bed"
python listRefExon.py refGene.txt.gz | sort -k1,1 -k2,2n -k3,3n - > refExon.bed

echo "bgzip -f refGene.bed"
bgzip -f refGene.bed

echo "bgzip -f refExon.bed"
bgzip -f refExon.bed

echo "tabix -f -p bed refGene.bed.gz"
tabix -f -p bed refGene.bed.gz 

echo "tabix -f -p bed refExon.bed.gz"
tabix -f -p bed refExon.bed.gz 


