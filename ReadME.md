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

1. Parsing breakpoint-containing and improperly aligned read pairs

```
GenomonSV parse sample.yaml param.yaml
```

2. Merging non-matched control panel breakpoint-containing read pairs
(for later filtering).

```
GenomonSV merge control.yaml mergedControl.bedpe.gz param.yaml                                        
```

3. Filtering and annotating candidate somatic structural variations

```
GenomonSV filt sample.yaml param.yaml
```

## Results


* **chr1**: chromosome for the 1st breakpoint
* **pos1**: coordinate for the 1st breakpoint
* **dir1**: direction of the 1st breakpoint
* **chr2**: chromosome for the 2nd breakpoint
* **pos2**: coordinate for the 2nd breakpoint
* **dir2**: direction of the 2nd breakpoint
* **inserted_seq**: inserted nucleotides within the breakpoints
* **variant_type**: type of the structural variation
* **gene1**: gene overlapping the 1st breakpoint
* **gene2**: gene overlapping the 2nd breakpoint
* **exon1**: exon overlapping the 1st breakpoint
* **exon2**: exon overlapping the 2nd breakpoint
* **#tumor_ref_read_pair**: #read_pairs not supporting the variant (reference read pairs) for the tumor sample
* **#tumor_var_read_pair**: #read_pairs supporting the variant (variant read paris) for the tumor sample
* **tumor_vaf**: frequency of variant read pairs for the tumor sample 
* **#control_ref_read_pair**: #read_pairs not supporting the variant for the matched control sample
* **#control_var_read_pair**: #read_pairs supporting the variant for the matched control sample
* **control_vaf**: frequency of variant read pairs for the matched control sample 
* **-log(Fisher P-value)**: minus common logarithm of p-value for the Fisher's exact text (on contingency table of (tumor v.s. matched control) and (reference v.s. variant read pairs)
* **non_matched_control_sample_with_max_junction**: sample name with the maximum number of junction read pairs
* **#max_non-matched_control_junction**: the maximum number of junction read pairs among non-matched control samples
* **max_over_hang1**: maximum overlang size of supporting read pairs from the 1st breakpoint
* **max_over_hang2**: maximum overlang size of supporting read pairs from the 2nd breakpoint

 
