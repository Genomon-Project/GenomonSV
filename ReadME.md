# Genomon SV

## Introduction

Genomon SV is a software for detecting somatic structural variations from cancer genome sequencing data.
Several characteristics of Genomon SV includes but not limited to;

1. Use both breakpoint-containing junction read pairs and improperly aligned read pairs for maximizing sensitivity
2. Various types of filters (e.g., use of non-matched normal control panels) are implemented for higher accuracy
3. Detection of short tandem duplications and mid-range deletions (10bp ~ 300bp) as well as larger structural variations such as translocations

## Dependency

### Python
Python (>= 2.7), `pysam (>= 0.8.1)`, `numpy`, `scipy`, `pyyaml` packages

### Software
tabix, bgzip, blat

## Install

```
git clone https://github.com/friend1ws/genomonSV.git
cd genomonSV
python setup.py build
python setup.py install
```
## Preparation

First, Genomon SV accept bam file aligned by `bwa mem` with -T0 option.
We do not guarantee the results for other cases.

Then, two types of configuration files (in yaml format) should be prepared.

1. sample.yaml
2. param.yaml
3. control.yaml

See sample files for description of each parameters.

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
* **#Tumor_Ref_Read_Pair**: #read_pairs not supporting the variant (reference read pairs) for the tumor sample
* **#Tumor_Var_Read_Pair**: #read_pairs supporting the variant (variant read paris) for the tumor sample
* **Tumor_VAF**: frequency of variant read pairs for the tumor sample 
* **#Control_Ref_Read_Pair**: #read_pairs not supporting the variant for the matched control sample
* **#Control_Var_Read_Pair**: #read_pairs supporting the variant for the matched control sample
* **Control_VAF**: frequency of variant read pairs for the matched control sample 
* **Minus_Log_Fisher_P_value**: minus common logarithm of p-value for the Fisher's exact text (on contingency table of (tumor v.s. matched control) and (reference v.s. variant read pairs)
* **Non-Matched_Control_Sample_With_Max_Junction**: sample name with the maximum number of junction read pairs
* **#Max_Non-Matched_Control_Junction**: the maximum number of junction read pairs among non-matched control samples
* **Max_Over_Hang_1**: maximum overlang size of supporting read pairs from the 1st breakpoint
* **Max_Over_Hang_2**: maximum overlang size of supporting read pairs from the 2nd breakpoint

 
