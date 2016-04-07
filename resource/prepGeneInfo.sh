#! /bin/bash

rm -rf refGene.txt.gz

wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz

echo "python listRefGene.py refGene.txt.gz | sort -k1,1 -k2,2n -k3,3n - > gene.bed"
python listRefGene.py refGene.txt.gz | sort -k1,1 -k2,2n -k3,3n - > gene.bed

echo "python listRefExon.py refGene.txt.gz | sort -k1,1 -k2,2n -k3,3n - > exon.bed"
python listRefExon.py refGene.txt.gz | sort -k1,1 -k2,2n -k3,3n - > exon.bed

echo "bgzip -f gene.bed"
bgzip -f gene.bed

echo "bgzip -f exon.bed"
bgzip -f exon.bed

echo "tabix -f -p bed gene.bed.gz"
tabix -f -p bed gene.bed.gz 

echo "tabix -f -p bed exon.bed.gz"
tabix -f -p bed exon.bed.gz 

