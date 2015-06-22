# Genomon SV

## Introduction

This is a repoistory for a new structural variation detection program, Genomon SV.

The name ebStruct is still tentative. But I'm considering to use some empirical bayesian framework
(Like EBCall, we perform some filter with Beta-Binomial artifact distribution from the many non-matched normal samples).


## Commands

## genomonSV_parse.sh

List up candidate SVs and their breakpoints (with improper read pairs).

## genomonSV_filt.sh

Perform various filters to reduce false positives.


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

 
