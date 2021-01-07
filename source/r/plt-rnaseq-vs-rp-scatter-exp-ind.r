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
# Program:  plt-rnaseq-vs-rp-scatter-exp.r 
# Version:  RPREP 1.0.0
# Author:   Travis L. Jensen, Johannes B. Goll, Sami R. Cherikh, and William F. Hooper
# Purpose:  Plot scatterplot of rnaseq and rp log2 expression
# Input:    
# Output:   N/A
#############################################################################################################

source('../r/init-analysis.r');

## specify directories
dir.lcpm = paste(res.dir,'lcpm',sep='/');

## align metadata
rp.mta = read.csv(paste0(dta.dir,'/sample_metadata.csv'),header=T,sep=',',stringsAsFactors=F)
rp.mta$ID = paste(rp.mta$subid,rp.mta$time,rp.mta$spct,rp.mta$trt,sep='_')
rs.mta = read.csv(paste0(gsub('ribosomal_profiling','rna_seq',dta.dir),'/sample_metadata.csv'),header=T,sep=',',stringsAsFactors=F)
rs.mta$ID = paste(rs.mta$subid,rs.mta$time,rs.mta$spct,rs.mta$trt,sep='_')

# limit and order to overlapping samples
rp.mta = rp.mta[rp.mta$ID %in% rs.mta$ID,]
rs.mta = rs.mta[match(rp.mta$ID,rs.mta$ID),]

if(nrow(rs.mta)>0) {
	
	########################
	#
	# Log2 expression scatter (Filtered)
	#
	#########################
	
	## import log2 expression data
	rp.lcpm.filt = read.csv(paste0(dir.lcpm,'/huh_alltp_lcpm_tmm_normalized_filtered.tab.gz'),header=T,sep='\t',stringsAsFactors=F)
	rs.lcpm.filt = read.csv(paste0(gsub('ribosomal_profiling','rna_seq',dir.lcpm),'/huh_alltp_lcpm_tmm_normalized_filtered.tab.gz'),header=T,sep='\t',stringsAsFactors=F)
	
	## plotting parameters
	par(mfrow=c(4,3),mar=c(3.5,3.8,2.5,0.2))
	
	## for each specimen type do
	for (j in 1:length(unique(rs.mta$spct))) {
		spc.flag = unique(rs.mta$spct)[j]
		spc.labl = unique(rs.mta$spctl)[j]
		
		## for ER or Cytosol do
		for (k in 1:length(unique(rs.mta$trt))) {
			trt.flag = unique(rs.mta$trt)[k]
			trt.labl = unique(rs.mta$trtl)[k]
		
			## for each timepoint do
			for (i in 1:length(unique(rs.mta$time))) {
				time.flag = sort(unique(rs.mta$time))[i]
				time.labl = unique(rs.mta$timel)[order(unique(rs.mta$time))][i]
				
				subs = unique(rp.mta[rp.mta$trt==trt.flag & rp.mta$spct==spc.flag & rp.mta$time==time.flag,'subid'])
				
				## for each subject do
				for (z in 1:length(subs)) {
					sub = subs[z]
					
					## subset metadata
					rp.mta.sel = rp.mta[rp.mta$trt==trt.flag & rp.mta$spct==spc.flag & rp.mta$time==time.flag & rp.mta$subid==sub,]
					rs.mta.sel = rs.mta[rs.mta$trt==trt.flag & rs.mta$spct==spc.flag & rs.mta$time==time.flag & rs.mta$subid==sub,]
					
					## intercecting features
					feat = intersect(rownames(rp.lcpm.filt),rownames(rs.lcpm.filt))
					
					## subset to sample
					if (nrow(rp.mta.sel)==1) {
						x = rp.lcpm.filt[match(feat,rownames(rp.lcpm.filt)),rp.mta.sel$samid]
						y = rs.lcpm.filt[match(feat,rownames(rs.lcpm.filt)),rs.mta.sel$samid]
					} else if (nrow(rp.mta.sel)>1) {
						x = rowMeans(rp.lcpm.filt[match(feat,rownames(rp.lcpm.filt)),rp.mta.sel$samid])
						y = rowMeans(rs.lcpm.filt[match(feat,rownames(rs.lcpm.filt)),rs.mta.sel$samid])
					} else {
						next;
					}
					
					## plot
					plot(x,y,main='',ylab='',xlab='',cex=0.5,col='#00000005',pch=15)
					mtext(side=1,line=2.6,expression('Mean Ribosomal Profiling Log'[2]*' CPM'),cex=0.85)
					mtext(side=2,line=2.1,expression('Mean RNA-Seq Log'[2]*' CPM'),cex=0.85)
					mtext(side=3,line=0.5,paste0(spc.labl,', ',trt.labl,', ',time.labl,', ',sub))
					
					## Add trendlines
					abline(1,1, col='grey50', lwd=1) 
					abline(lm(y~x), col='orange', lwd=1.4) 
					lines(lowess(y~x), col='forestgreen', lwd=1.4) 	
					
					## Add trendline legend
					legend('bottomright',legend=c('Linear','LOWESS','1:1'),col=c('orange','forestgreen','grey50'), lwd=1.2,bty='n',lty=c(1,1))
					
					## Add legend with linear/r/rs/n
					r  = cor.test(x,y,method='pearson') 
					rs = cor(x,y,method='spearman')
					lm.coef = round(lm(y~x)$coefficients,2)
					legend.text = c(paste('y = ',lm.coef[2],'x ',if(lm.coef[1]>=0) {'+ '},lm.coef[1],sep=''),
							paste0('r = ',round(r$estimate,3)),
							paste0('rs = ',round(rs,3)),
							paste0('genes = ',length(x)))
					legend('topleft',bty='n',legend=legend.text)
				}
			}
		}
	}
}