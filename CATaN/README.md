# CATaN

<!-- badges: start -->
<!-- badges: end -->

**C**anonical Correlation **A**nalysis of **T**ranscriptome **a**nd
TF-gene regulatory **N**etworks

## Overview

CATaN performs Canonical Correlation Analysis (CCA) between transcription
factor (TF) gene-regulatory networks matrix (TF-GRN) and transcriptome data to identify shared
regulatory axes (Canonical Components). CC scores are then aligned to genomic
SNP positions for use as annotations in stratified LD Score Regression
(S-LDSC).

## Installation

```r
# Install development version from GitHub
install.packages("BiocManager")
BiocManager::install("htakahashi1/CATaN")
```

## Quick start

```r
library(CATaN)

# Step 1: Filter and normalise
filtered <- filter_low_expression(counts)
mats <- prepare_matrices(filtered, tf_matrix)

# Step 2: Run CCA
result <- run_cca(mats$tf, mats$expr, n_cc = 10L)

# Step 3: Align to SNPs and export for S-LDSC
annot <- align_cc_to_snps(result, peak_dir, snp_gr)
annot <- extract_top_bottom_snps(annot)
export_for_sldsc(annot, "output/")
```

See `vignette("CATaN_tutorial")` for a full walkthrough.

## Citation

If you use CATaN in your research, please cite:

> <Preparing>
