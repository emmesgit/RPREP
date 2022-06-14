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
# Program:  read-length-distribution-figures.r 
# Version:  RPREP 1.1.0
# Author:   Travis L. Jensen, Johannes B. Goll, Sami R. Cherikh, and William F. Hooper
# Purpose:  Generate density plot using read length data calculated in RP preprocessing, 
#           colored by pbmc/rna amount. Save/load density object to avoid recomputing it 
#           every time
# Input:    data/length_distribution/<samid>_<stage>_counts.RData OR
#           data/length_distribution/<samid>_<stage>_density.RData
# Output:   N/A
#############################################################################################################

source('../r/init-analysis.r');
options(scipen=999)

## Input directory
in.dir = paste(dta.dir, '/rlength_distribution', sep='')

## Plot pre & post read length distributions
lengthFilterFlags = c('pre', 'post')
lengthFilterLabls = c('Before read length filtering', 'After read length filtering')

## plotting parameters
par(mfrow=c(2,1),oma=c(0,0,0,0),mar=c(3,3,2,0.5))
    
## For each filter flag
for (f in 1:length(lengthFilterFlags)){
    lengthFilterFlag = lengthFilterFlags[f]
    lengthFilterLabl = paste0(lengthFilterLabls[f])
    
    ## For each sample
    for (s in 1:length(mta$samid)) {
        samid = mta$samid[s]
        
        in.file.pre = paste0(in.dir,'/',samid,'_1_',lengthFilterFlag,'_counts.RData')
        if (!file.exists(in.file.pre)) next
          
        load(in.file.pre)
        
        ## Only plot panel & label axes once 
        if (s == 1) {
            plot(-1,-1,
                 type='l',
                 ylim=c(0,0.5),
                 xlim=c(0,100),
                 xlab="", 
                 cex.main=0.8,
                 ylab="",
                 lwd=0.8,
                 main=lengthFilterLabl, cex.axis=0.7)
            mtext("Length (nt)",side=1,line=2,cex=0.7)
            mtext("Frequency",side=2,line=2,cex=0.7)
        }
        
        ## Plot lines with slight transparency
        lines(d, col=makeTransparent(mta$spctc[s], alpha=80), lwd=1.5)
    }
}
