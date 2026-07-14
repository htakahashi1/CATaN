# CATaN

<!-- badges: start -->
<!-- badges: end -->

**C**anonical Correlation **A**nalysis of **T**ranscriptome **a**nd
TF-gene regulatory **N**etworks

## Overview

CATaN performs Canonical Correlation Analysis (CCA) between transcription
factor (TF) gene-regulatory networks matrix (TF-GRN) and transcriptome data to identify shared
regulatory axes (Canonical Components). CC scores are then aligned to genomic
SNP positions to make binary annotations for stratified LD Score Regression
(S-LDSC).

## Installation

```r
# Install development version from GitHub
install.packages("BiocManager")
BiocManager::install("htakahashi1/CATaN")
```

## Tutorial

See https://htakahashi1.github.io/CATaN/ for a full walkthrough.

## Citation

If you use CATaN in your research, please cite:

> <Preparing>
