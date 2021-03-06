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
# Program:  run-report.sh
# Version:  RPREP 1.0.1
# Author:   Travis L. Jensen, Johannes B. Goll, Sami R. Cherikh, and William F. Hooper
# Purpose:  script to run reporting only steps (knitR/latex based results reporting)
# Input:    N/A
# Output:   N/A
#############################################################################################################

## Exit script upon the failure of any child scripts
#set -e 
#set -o pipefail

ACD=$1
SRCDIR=$2/source

## construct report
rm $ACD/report/cache/* $ACD/report/figs/* $ACD/report/tables/* || true
Rscript -e "library(knitr); setwd('$ACD/report'); knit('$SRCDIR/knitr/rprep-report.Rnw')";
cd $ACD/report;
pdflatex -interaction=batchmode rprep-report.tex;
pdflatex -interaction=batchmode rprep-report.tex;

## convert PDF to PNG
cd $ACD/report/figs
for f in *.pdf; do
    gs -q -dBATCH -dNOPAUSE -sDEVICE=pngalpha -r300x300 -sOutputFile=${f%.pdf}.png $f
done

