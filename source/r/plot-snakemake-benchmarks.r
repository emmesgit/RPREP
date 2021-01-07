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
# Program:  plot-snakemake-benchmarks.r 
# Version:  RPREP 1.0.0
# Author:   Travis L. Jensen, Johannes B. Goll, Sami R. Cherikh, and William F. Hooper
# Purpose:  Plot snakemake benchmarks
# Input:    argv[1] : Path to merged snakemake benchmarks
# Output:   argv[1]/benchmark_results.pdf
#############################################################################################################

libs = c('RColorBrewer')
suppressPackageStartupMessages(lapply(libs, require, character.only=T))

## Init plotting units/titles 
units  = c('Seconds', 'MB','MB','MB','MB','MB','MB','%')
titles = c('Runtime', 'Max Resident Set Size', 'Max Virtual Memory Size', 
           'Max Unique Set Size','Max Proportional Set Size', 'I/O In',
           'I/O Out', 'Mean Load')

## Set input path
argv = commandArgs(trailingOnly=T)[1]
if (!dir.exists(argv)) stop('Directory',argv,'does not exist!')
setwd(argv)


## Get list of files and merge
files = list.files(argv, pattern='.tab$')
rules = sapply(files, function(x) unlist(strsplit(x,'.',fixed=T))[1])
res = c()
for (f in 1:length(files)) {
    dta      = read.csv(files[f], h=T, stringsAsFactors=F, sep='\t')
    dta$rule = rules[f]
    res      = rbind(res, dta)
}

## Remove h:m:s field
res = as.data.frame(res[,-2])

## If there are more than 4 rules, don't try to squeeze two plots per row
if (length(files) < 5){
    par(mfrow=c(4,2), mar=c(2,4,4,2), mgp=c(2,0.75,0))
} else {
    par(mfrow=c(4,1), mar=c(2,4,4,2), mgp=c(2,0.75,0))
}

## Plot each metric
for (i in 1:(ncol(res)-1)){
    col = colnames(res)[i]

    boxplot(tapply(as.numeric(res[[col]]), res$rule, c),
            col=terrain.colors(length(files)), 
            main=titles[i],
            cex.axis=0.5,
            xlab='Rule',
            ylab=units[i])
}


