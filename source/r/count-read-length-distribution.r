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
# Program:  count-read-length-distribution.r
# Version:  RPREP 1.0.1
# Author:   Travis L. Jensen, Johannes B. Goll, Sami R. Cherikh, and William F. Hooper
# Purpose:  Given a sequence of integers corresponding to read length, plot size distribution
# Input:    (stdin from awk 1-liner) sequence of integers corresponding to read length in a fasta file
# Output:   pre.dir/stage1/rlength_distribution/<samid>_<pre|post|raw>_counts.RData
#           pre.dir/stage1/rlength_distribution/<samid>_<pre|post|raw>_counts.RData
#############################################################################################################

rm(list=ls(all=TRUE))

## Grab command line arguments, make sure all 4 are present:
## argv[1] : sample_id (used to name file output)
## argv[2] : pdf output path 
## argv[3] : workspace image output path 
## argv[4] : raw data, pre-, or post- length-filtered data
argv = commandArgs(trailingOnly=TRUE)

if (length(argv) != 4){
    ## Print error message and exit with non-zero value
    print('Error: One or more command line arguments are missing. Script needs samid, pdf output path, RData output path, data type')
    q("no", 1, FALSE)
}

if (argv[4] == 'pre') {
    order = "Prefilter"
} else if (argv[4] == 'post') {
    order = "Postfilter"
} else if (argv[4] == 'raw') {
    order = "Unprocessed"
}

sample_id = argv[1]
pdf_out = argv[2]
data_out = argv[3]



## Grab counts from stdin 
counts = read.table(file("stdin"), h=F, stringsAsFactors=F)

## Plot histogram to pdf
pdf(file=pdf_out)
main.title = paste(order, "read length distribution for", sample_id, sep=" ")
if (order == "Postfilter"){
    d = density(counts$V1, adjust=10)
    hist(counts$V1, plot=T, main=main.title, xlab="Read length (nt)", freq=FALSE, breaks=(20:40), col='grey', ylim=c(0,max(d$y)))
} else {
    d = density(counts$V1, adjust=10)
    hist(counts$V1, plot=T, main=main.title, xlab="Read length (nt)", freq=FALSE, col='grey', ylim=c(0,max(d$y)))
}
lines(d, col="blue")
dev.off()

## Delete counts for faster loading later, smaller filesize
rm(counts)

## Write out Rdata file
save.image(file=data_out, compress=TRUE)