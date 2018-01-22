# Genomon SV

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Build Status](https://travis-ci.org/Genomon-Project/GenomonSV.svg?branch=devel)](https://travis-ci.org/Genomon-Project/GenomonSV)

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
pip install genomon-sv
```
For the last command, you may need to add --user if using a shared computing cluster.

## Preparation

Before using GenomonSV, 3 preparation steps are required

### Install required softwares and add them to the PATH

GenomonSV uses tabix, bgzip (which are part of HTSlib projects) and blat inside the programs, 
assuming those are installed and the pathes are already added to the running environment.
Please install and add them to the PATH.


### Prepare bam files

Genomon SV accept just bam file aligned by `bwa mem` with -T0 option.
We do not guarantee the results for other cases.
Also, we assume that the sequencing data is paired-end. All the single-end reads are ignored in the program.


## Commands

### parse

Parsing breakpoint-containing and improperly aligned read pairs
```
GenomonSV parse [-h] [--debug]
                [--junction_abnormal_insert_size JUNCTION_ABNORMAL_INSERT_SIZE]
                [--junction_min_major_clipping_size JUNCTION_MIN_MAJOR_CLIPPING_SIZE]
                [--junction_max_minor_clipping_size JUNCTION_MAX_MINOR_CLIPPING_SIZE]
                [--junction_check_margin_size JUNCTION_CHECK_MARGIN_SIZE]
                [--improper_abnormal_insert_size IMPROPER_ABNORMAL_INSERT_SIZE]
                [--improper_min_mapping_qual IMPROPER_MIN_MAPPING_QUAL]
                [--improper_max_clipping_size IMPROPER_MAX_CLIPPING_SIZE]
                [--junction_dist_margin JUNCTION_DIST_MARGIN]
                [--junction_opposite_dist_margin_margin JUNCTION_OPPOSITE_DIST_MARGIN_MARGIN]
                [--improper_check_margin_size IMPROPER_CHECK_MARGIN_SIZE]
                input.bam output_prefix
```
- **input.bam**: Path to input indexed bam file
- **output_prefix**: Output file prefix
See the help (``GenomonSV parse -h``) for other options. 
But I believe the default option is 
enough for typical illumina sequence data (with the length of 50bp ~ 200bp).

After successful completion, you will find possible breakpoint regions evidenced by
breakpoint containing read pairs ({output_prefix}/junction.clustered.bedpe.gz(.tbi))
and improperly aligned read pairs ({output_prefix}/improper.clustered.bedpe.gz(.tbi))).


### merge

Merging non-matched control panel breakpoint-containing read pairs.
This step picks up germline and artifacts breakpoints (e.g., black-list breakpoints) for later filtering steps,
typically using several control samples.
We strongly believe this step is crucial for improving accuracy of somatic structural variation calling.

```
GenomonSV merge [-h] [--debug]
                [--merge_check_margin_size MERGE_CHECK_MARGIN_SIZE]
                control_info.txt merge_output_file                                     
```
- **control_info.txt**: Tab-delimited file on non-matched control. 
The 1st column is sample label for each breakpoint information file, and can be freely specified.
The 2nd column is the output_prefix generated at the above parse stage
(GenomonSV merge program assumes each {output_prefix}.junction.clustered.bedpe.gz file is already generated).
- **merge_output_file**: Output merged breakpoint information file

### filt
Filtering and annotating candidate somatic structural variation.

```
GenomonSV filt [-h] [--matched_control_bam matched_control.bam]
                      [--non_matched_control_junction merged.junction.control.bedpe.gz]
                      [--matched_control_label MATCHED_CONTROL_LABEL]
                      [--genome_id {hg19,hg38,mm10}] [--grc] [--debug]
                      [--min_junc_num MIN_JUNC_NUM]
                      [--min_sv_size MIN_SV_SIZE]
                      [--min_inversion_size MIN_INVERSION_SIZE]
                      [--control_panel_num_thres CONTROL_PANEL_NUM_THRES]
                      [--control_panel_check_margin CONTROL_PANEL_CHECK_MARGIN]
                      [--min_support_num MIN_SUPPORT_NUM]
                      [--min_mapping_qual MIN_MAPPING_QUAL]
                      [--min_overhang_size MIN_OVERHANG_SIZE]
                      [--close_check_margin CLOSE_CHECK_MARGIN]
                      [--close_check_thres CLOSE_CHECK_THRES]
                      [--max_depth MAX_DEPTH] [--search_length SEARCH_LENGTH]
                      [--search_margin SEARCH_MARGIN]
                      [--split_refernece_thres SPLIT_REFERNECE_THRES]
                      [--validate_sequence_length VALIDATE_SEQUENCE_LENGTH]
                      [--short_tandem_reapeat_thres SHORT_TANDEM_REAPEAT_THRES]
                      [--blat_option BLAT_OPTION]
                      [--min_tumor_variant_read_pair MIN_TUMOR_VARIANT_READ_PAIR]
                      [--min_tumor_allele_freq MIN_TUMOR_ALLELE_FREQ]
                      [--max_control_variant_read_pair MAX_CONTROL_VARIANT_READ_PAIR]
                      [--max_control_allele_freq MAX_CONTROL_ALLELE_FREQ]
                      [--max_fisher_pvalue MAX_FISHER_PVALUE]
                      input.bam output_prefix reference.fa
```

- **input.bam**: Path to input indexed bam file
- **output_prefix**: Output file prefix (assuming files generated at the parse step already exist).
- **reference.fa**: Path to reference genome sequence (fasta format) used for alignment.
- **annotation_dir**: Path to the gene and exon information file (in which gene.bed.gz(.tbi) and exon.bed.gz(.tbi) exist).

The following options are not mandatory, but we strongly believe is necessary for improved results.
- **--matched_control_bam**: The path to matched control bam file. when this is specified, GenomonSV performs Fisher's exact test on the numbers of variant and non-variant read pairs between tumor and control and remove those with high p-value. we believe this filtering step is crucial for removing massive false positives and highly recommend to use this. even when matched control sequence data is not available, using *dummy* mathced control may be helpful.
- **--non_matched_control_junction**: The path to merged junction files created at the merge step. for each structural variation candidate, when the number of supporing junction read pairs shared by any of the pooled control samples is equal to or more than the specified threshould (**control_panel_num_thres** option), then that candidate is filtered out. 
- **--matched_control_label**: In the above filtering step using non-matched control junctions, ignore the sample speficied by this option. Typically, matched control sample label is specified for avoiding filtering true positives because of tumor cell contamination in the control sample.

See the help (``GenomonSV filt -h``) for other options.
You may want to tune up **min_junc_num**, **min_support_num**, **min_overhang_size**, **max_depth**, **min_tumor_variant_read_pair**,
**min_tumor_allele_freq**, **max_control_variant_read_pair**, **max_control_allele_freq**, **-max_fisher_pvalue**
depending on your sequencing depth and tumor purity.

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

 
