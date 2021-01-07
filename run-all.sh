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
# Program:  run-all.sh
# Version:  RPREP 1.0.0
# Author:   Travis L. Jensen, Johannes B. Goll, Sami R. Cherikh, and William F. Hooper
# Purpose:  run ll steps of RPREP
# Input:    N/A
# Output:   N/A
#############################################################################################################

## Locate Source Directory from root
SRCDIR="$(cd `dirname $0` && pwd)/source";

## dont forget to specify minimally "--config XXX", "--threads XXX", and "--log XXX" with this shell script command:
## for example: sh run-all.sh --config case-study/config-dengue.xlsx --threads 32 --log run-all.log
python3 $SRCDIR/python/rprep --assay rna_seq --run-preprocessing --gene-body-coverage $@
python3 $SRCDIR/python/rprep --assay ribosomal_profiling --run-preprocessing --gene-body-coverage $@
python3 $SRCDIR/python/rprep --assay rna_seq --run-analysis $@
python3 $SRCDIR/python/rprep --assay ribosomal_profiling --run-analysis $@
python3 $SRCDIR/python/rprep --report $@
