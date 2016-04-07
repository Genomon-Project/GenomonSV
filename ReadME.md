# Genomon SV

## Introduction

Genomon SV is a software for detecting somatic structural variations from cancer genome sequencing data.
Several characteristics of Genomon SV includes but not limited to;

1. Use both breakpoint-containing junction read pairs and improperly aligned read pairs for maximizing sensitivity
2. Various types of filters (e.g., use of non-matched normal control panels) are implemented for higher accuracy
3. Detection of short tandem duplications and mid-range deletions (10bp ~ 300bp) as well as larger structural variations such as translocations

## Dependency

### Python
Python (>= 2.7), `pysam (>= 0.8.1)`, `numpy`, `scipy` packages

### Software
[hstlib](http://www.htslib.org), [blat](http://hgdownload.cse.ucsc.edu/admin/exe/)

## Install

```
git clone https://github.com/friend1ws/GeomonSV.git
cd GenomonSV
python setup.py build install
```
For the last command, you may need to add --user if using a shared computing cluster.

## Preparation

Before using GenomonSV, 3 preparation steps are required

### Install required softwares and add them to the PATH

GenomonSV uses tabix, bgzip (which are part of HTSlib projects) and blat inside the programs, 
assuming those are installed and the pathes are already added to the running environment.
Please install and add them to the PATH.

### Prepare annotation files

At the last step, GenomonSV add annotation information to each structural variation candidates.
We need bgzip compressed bed format gene and exon information files named as gene.bed.gz and exon.bed.gz,
with tabix index files (gene.bed.gz.tbi and exon.bed.gz.tbi).

We prepared a sample program for preparing them

```
cd GenomonSV/resource
bash prepGeneInfo.sh
```

### Prepare bam files

Genomon SV accept just bam file aligned by `bwa mem` with -T0 option.
We do not guarantee the results for other cases.


## Commands

* Parsing breakpoint-containing and improperly aligned read pairs

```
GenomonSV parse sample.yaml param.yaml
```

* Merging non-matched control panel breakpoint-containing read pairs
(for later filtering).

```
GenomonSV merge control.yaml mergedControl.bedpe.gz param.yaml                                        
```

* Filtering and annotating candidate somatic structural variations

```
GenomonSV filt sample.yaml param.yaml
```

## Results


* **Chr_1**: chromosome for the 1st breakpoint
* **Pos_1**: coordinate for the 1st breakpoint
* **Dir_1**: direction of the 1st breakpoint
* **Chr_2**: chromosome for the 2nd breakpoint
* **Pos_2**: coordinate for the 2nd breakpoint
* **Dir_2**: direction of the 2nd breakpoint
* **Inserted_Seq**: inserted nucleotides within the breakpoints
* **Variant_Type**: type of the structural variation
* **Gene_1**: gene overlapping the 1st breakpoint
* **Gene_2**: gene overlapping the 2nd breakpoint
* **Exon_1**: exon overlapping the 1st breakpoint
* **Exon_2**: exon overlapping the 2nd breakpoint
* **Num_Tumor_Ref_Read_Pair**: #read_pairs not supporting the variant (reference read pairs) for the tumor sample
* **Num_Tumor_Var_Read_Pair**: #read_pairs supporting the variant (variant read paris) for the tumor sample
* **Tumor_VAF**: frequency of variant read pairs for the tumor sample 
* **Num_Control_Ref_Read_Pair**: #read_pairs not supporting the variant for the matched control sample
* **Num_Control_Var_Read_Pair**: #read_pairs supporting the variant for the matched control sample
* **Control_VAF**: frequency of variant read pairs for the matched control sample 
* **Minus_Log_Fisher_P_value**: minus common logarithm of p-value for the Fisher's exact text (on contingency table of (tumor v.s. matched control) and (reference v.s. variant read pairs)
* **Non-Matched_Control_Sample_With_Max_Junction**: sample name with the maximum number of junction read pairs
* **Num_Max_Non-Matched_Control_Junction**: the maximum number of junction read pairs among non-matched control samples
* **Max_Over_Hang_1**: maximum overlang size of supporting read pairs from the 1st breakpoint
* **Max_Over_Hang_2**: maximum overlang size of supporting read pairs from the 2nd breakpoint

 
