#!/bin/bash
#############################################################################################################
# RPREP: Ribosomal Profiling Reports, an open-source cloud-enabled framework for reproducible
# Ribosomal Profiling data processing, analysis, and result reporting
# 
# https://github.com/emmesgit/RPREP
#
# Copyright (C) 2020 The Emmes Company L.L.C. 
# 
# This program is free software that contains third party software subject to various licenses, 
# namely, the GNU General Public License version 3 (or later), the GNU Affero General Public License 
# version 3 (or later), and the LaTeX Project Public License v.1.3(c). A list of the software contained 
# in this program, including the applicable licenses, can be accessed here: 
# 
# https://github.com/emmesgit/RPREP/blob/master/SOFTWARE.xlsx  
# 
# You can redistribute and/or modify this program, including its components, only under the terms of 
# the applicable license(s).  
#
# To cite this software, please reference <DOI>
#
# This program is distributed in the hope that it will be useful, but "as is," WITHOUT ANY WARRANTY; 
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
#
# Program:  run-analysis.sh
# Version:  RPREP 1.0.0
# Author:   Travis L. Jensen, Johannes B. Goll, Sami R. Cherikh, and William F. Hooper
# Purpose:  script to run analysis only steps (R based normalization, 
#                PCA and MDS, DE gene Identification, gene clusters, and GSEA)
# Input:    N/A
# Output:   N/A
#############################################################################################################

## Exit script upon the failure of any child scripts
set -euo pipefail

## Locate Source Directory from root
SRCDIR="$(cd `dirname $0` && pwd)/source";

## get preprocessing directory, analysis directory, workflow configuration and metadata csv file locations
WCD=`head -1 $SRCDIR/../dir.csv`;
ACD=`head -2 $SRCDIR/../dir.csv | tail -1`;
ASSAY=`tail -2 $SRCDIR/../dir.csv | head -1`;

## migrate to R directory
cd $SRCDIR/r

## update sample metadata -- add bam statistics
echo 'Adding RSeQC metrics to sample metadata'
Rscript 00-qc-normalization/init-sample-meta-data.r $ACD/data/$ASSAY $ACD/data/$ASSAY $ACD/data/$ASSAY 

## import annotations from ensembl
echo 'Downloading annotations from Ensembl'
Rscript 00-qc-normalization/init-export-annotations.r 

## perform TMM normalization
echo 'Running TMM normalization'
Rscript 00-qc-normalization/init-tmm-normalization-fragments.r 

## calculate log fold changes
echo 'Computing log fold changes'
Rscript 00-qc-normalization/init-log-cpm-fold-change-from-baseline.r 

## calculate PCA, Euclidean, and Spearman distances (Run in parallel)
echo 'Running PCA; Computing Euclidean and Spearman distances'
Rscript 01-bias-confounding-effects/init-euclid-dist.r &
Rscript 01-bias-confounding-effects/init-spearman-dist.r &
Rscript 01-bias-confounding-effects/init-pca.r 

## discover significantly differentially expressed genes (SDEG)
echo 'Running differential gene analysis'
Rscript 02-sdeg-identification/init-edgeR-glm-model.r 

## initialize pvclusters
echo 'Clustering differentially expressed genes'
Rscript 03-sdeg-clusters/init-pvclusters.r 

## run gene set enrichment analysis (GSEA)
echo 'Running pathway enrichment analysis'
Rscript 04-sdeg-organization-known-modules/init-gsea-sampling.r;

