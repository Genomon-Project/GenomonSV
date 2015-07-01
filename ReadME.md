# Genomon SV

## Introduction

Genomon SV is a software for detecting somatic structural variations from cancer genome sequencing data.
Several characteristics of Genomon SV includes but not limited to;

1. Use both breakpoint-containing junction read pairs and improperly aligned read pairs for maximizing sensitivity
2. Various types of filters (e.g., use of non-matched normal control panels) are implemented for higher accuracy
3. Detection of short tandem duplications and mid-range deletions (10bp ~ 300bp) as well as larger structural variations such as translocations

## Dependency

Python (>= 2.7), `pysam (>= 0.8.1)`, `numpy`, `scipy`, `pyyaml`, `Biopython` packages

## Install

```
python setup.py build
python setup.py install
```
## Preparation

First, Genomon SV accept bam file aligned by `bwa mem` with -T0 option.
We do not guarantee the results for other cases.

Then, two types of configuration files (in yaml format) should be prepared.

1. sample.yaml
2. param.yaml

See sample files for description of each parameters.

## Commands

1. Parsing breakpoint-containing and improperly aligned read pairs

```
GenomonSV parse sample.yaml param.yaml
```

2. Filtering and annotating candidate somatic structural variations

```
GenomonSV filt sample.yaml param.yaml
```

## Results

1. chromosome for the 1st breakpoint
1. coordinate for the 1st breakpoint
1. direction of the 1st breakpoint
1. chromosome for the 2nd breakpoint
1. coordinate for the 2nd breakpoint
1. direction of the 2nd breakpoint
1. inserted nucleotides within the breakpoints
1. type of the structural variation
1. gene overlapping the 1st breakpoint
1. gene overlapping the 2nd breakpoint
1. exon overlapping the 1st breakpoint
1. exon overlapping the 2nd breakpoint
1. #read_pairs not supporting the variant (reference read pairs) for the tumor sample
1. #read_pairs supporting the variant (variant read paris) for the tumor sample
1. #read_pairs not supporting the variant for the matched control sample
1. #read_pairs supporting the variant for the matched control sample
1. frequency of variant read pairs for the tumor sample 
1. frequency of variant read pairs for the matched control sample 
1. p-value for the Fisher's exact text (on contingency table of (tumor v.s. matched control) and (reference v.s. variant read pairs)

 
