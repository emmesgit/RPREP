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
# Program:  gene-coverage-bias-plots.r 
# Version:  RPREP 1.0.1
# Author:   Travis L. Jensen, Johannes B. Goll, Sami R. Cherikh, and William F. Hooper
# Purpose:  Merge sample gene body coverage data and generate curve/heatmap
# Input:    data/gene_body_coverage/bias_plots/<sample>.geneBodyCoverage.txt
# Output:   N/A
#############################################################################################################

source('../r/init-analysis.r')
source('../r/heatmap4.r')
options(scipen=999)

## Set input dir
in.dir = paste0(dta.dir,'/bias_plots')

## Reorder metadata, set colors
mta = mta[order(mta$subid,mta$time),]

## Merge data, taking into account that samples could be missing 
dta = data.frame()
missing.samples = c()

for (s in 1:length(mta$samid)){
    in.file = paste(in.dir,'/',mta$samid[s],'.geneBodyCoverage.txt',sep='')
    
    if (file.exists(in.file)) {
        tmp = read.table(in.file, h=T, stringsAsFactors=F)[,2:101]
        dta = rbind(dta, tmp)    
    } else {
        missing.samples = c(missing.samples, s)        
    }
}



## Don't attempt to plot if no coverage data is available
if (!length(missing.samples) == nrow(mta)) {
    
    if (length(missing.samples) > 0) {
        rownames(dta) = mta$samid[-which(mta$samid %in% missing.samples)]    
    } else {
        rownames(dta) = mta$samid
    }
    
    ## Calculate percent coverage from # reads
    dta = dta / apply(dta,1,max)
    
    ## Plotting parameters
    par(mfrow=c(1,2),oma=c(0,0,0,0))
    rowCols  = mta$spctc[mta$samid %in% rownames(dta)]
    rowSpc   = mta$spctl[mta$samid %in% rownames(dta)]
    
    
    
    #########################
    #        Line plot        #
    #########################
    
    par(mar=c(3,4,3,1))
    plot(-1,-1,type='l',ylim=c(0,1),xlim=c(1,100),xlab="", cex.main=0.8,
            ylab="",lwd=0.8,main='Gene Body Coverage\n')
    mtext("Gene body percentile (5'->3')",side=1,line=2,cex=0.7)
    mtext("Coverage",side=2,line=2,cex=0.7)
    
    for (i in 1:nrow(dta)) {
        y = as.numeric((dta[i,]))
        lines(1:100,y,type='l',col=rowCols[i])
    }
    
    ## Add specimen types legend
    legend('bottom', legend=unique(rowSpc), fill=rowCols, bty='n', horiz=T, cex=0.7)
    
    
    
    #########################
    #         Heatmap        #
    #########################
    
    ## Need at least 2 samples to plot a heatmap
    if (nrow(dta) > 1) {
        axis.loc = c(1,10,20,30,40,50,60,70,80,90,100)
        
        par(mar=c(3,1,3,5))
        rc = cm.colors(100)
        heatmap4(as.matrix(dta), scale=c("none"),labRow = rownames(dta) ,Colv = NA,Rowv = NA,
                labCol=NA,col=cm.colors(256),ColSideColors = rc,cexRow=0.5, 
                cexCol=1,xlab='',
                add.expr=x_axis_expr <- axis(side=1,at=axis.loc, labels=as.character(axis.loc)))
        title(main='Gene Body Coverage\n',cex.main=0.8)
        mtext("Gene body percentile (5'->3')",side=1,line=2,cex=0.7)
    }
}
