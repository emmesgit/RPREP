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
# Program:  get-rrna-fastas.r 
# Version:  RPREP 1.1.0
# Author:   Travis L. Jensen, Johannes B. Goll, Sami R. Cherikh, and William F. Hooper
# Purpose:  Pulls human rRNA, mt-rRNA, mt-tRNA from Ensembl using biomaRt, writes out fasta
# Input:    argv[1] : output fasta location (usually pre.dir/stage1/index/ens_xrna.fasta)
#           argv[2] : Ensembl version
# Output:   Fasta specified in argv[1]
#           dir(argv[1])/ens-rrna-version.txt (logfile)
#############################################################################################################

## Reset variables, load library
rm(list=ls(all=TRUE))

libs = c('biomaRt','stringr')
lapply(libs, require, character.only=TRUE)

## argv[1] : FASTA output dir+filename. Logfile will be placed in the same dir 
## argv[2] : Ensembl version number
argv = commandArgs(trailingOnly=TRUE)

## Don't run unless both arguments are present
if (length(argv) !=2){
    ## Print error message and exit with non-zero value
    print('Error: One or more command line arguments are missing. Script needs FASTA output file & ensembl version specified')
    q("no", 1, FALSE)
}

outfile = argv[1]
ensembl.version = as.numeric(argv[2])
logfile = paste(dirname(outfile),'/ensembl-rrna-version.txt' ,sep='')

## Set mart as ensembl, dataset as Homo sapiens
## Grab ensembl gene id's for rRNA, mt-rRNA, mt-tRNA 
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version=ensembl.version)
att=c('ensembl_gene_id', 'gene_biotype')
fil = c('biotype')
val = c('rRNA', 'Mt_rRNA', 'Mt_tRNA', 'Mt_tRNA_pseudogene', 'tRNA_pseudogene', 'rRNA_pseudogene')
ensembl_ids = getBM(attributes=att, filters=fil, values=val, mart=ensembl)

## Grab sequences
sequence = getSequence(id=ensembl_ids$ensembl_gene_id, type="ensembl_gene_id", mart=ensembl, seqType='cdna')

## getSequence() doesn't always return genes in the same order.. sort by ensembl ID before exporting
sequence = sequence[order(sequence$ensembl_gene_id), ]

## If the file already exists, delete it.. no option to force overwrite for exportFasta(), it just appends
if (file.exists(outfile)) system(paste('rm',outfile))
exportFASTA(sequence, outfile)


## Write log with Ensembl version, date of download
sink(logfile)
cat(paste("Ensembl version ", ensembl.version, '\n', sep=''))
cat(paste("Date downloaded: ", date(), sep=''))
sink()

