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
# Program:  translation-efficiency-tables.r
# Version:  RPREP 1.0.1
# Author:   Travis L. Jensen, Johannes B. Goll, Sami R. Cherikh, and William F. Hooper
# Purpose:  Plot table identifying similarities/differences between RP/RNA-Seq SDEGs
# Input:    glm/<trt>_<day>_glm_sig.tab.gz
#           data/annot/gr38_filtered_gene_annotations.tab.gz
#           lcpm/alld_lcpm_tmm_normalized_filtered.tab.gz
# Output:   N/A
#############################################################################################################

source('../r/init-analysis.r')

## Dont print in scientific notation -- use 6 digits following decomal
options("scipen"=100, "digits"=6)

## Load gene annotations
ano.file = paste0(dta.dir,'/annot/filtered_gene_annotations.tab.gz')
gen      = read.csv(ano.file, stringsAsFactors=F, h=T ,sep='\t')



## Get significant genes & log fold changes for each assay, return in a list
getGenes = function(in.file, dir=bas.dir, assays=assayFlags) {
    res = vector('list', 2)
    
    ## Collect genes & fold changes for each assay
    for (a in 1:length(assays)) {
        assay = assays[a]
        f       = paste0(bas.dir,'/analysis/',assay,'/glm/',in.file)
        if (file.exists(f)) {
            x = read.csv(f, h=T, stringsAsFactors=F, sep='\t')
            res[[assay]] = x[,'logFC']
            names(res[[assay]]) = rownames(x)
        } else {
            res[[assay]] = c()
        }
    }

    return(res)
}



constructRegString = function(gen, lfc, sig.genes) {
    if (is.na(lfc)) {
        return('LE')
    } else if (!gen %in% sig.genes) {
        return('NS')
    } else {
        return(ifelse(lfc > 0, '+', '-'))
    }
}


## For each specimen type
for(s in 1:length(spcFlags)) {
    spcFlag = spcFlags[s];
    spcLabl = spcLabls[s];
    
    ## For each treatment group
    for (v in 1:length(trtFlags)) {
        trtFlag = trtFlags[v];
        trtLabl = trtLabls[v];
        
        ## For each post-baseline timepoint
        for(t in 1:length(postb.times)) {
            time = postb.times[t];
            timel = postb.times.l[t];
            
            ## Get significant genes for each assay
            in.file = paste0(spcFlag,'_',trtFlag,'_tp',time,'_glm_sig.tab.gz')
            gen.sig = getGenes(in.file=in.file)
            
            ## Get all genes
            in.file = paste0(spcFlag,'_',trtFlag,'_tp',time,'_glm_all.tab.gz')
            gen.all = getGenes(in.file=in.file)
            
            ## Compute union of significant genes across both assays
            gen.sig.union = c()
            for (g in 1:length(gen.sig)) gen.sig.union = union(gen.sig.union, names(gen.sig[[g]]))
            
            
            ## For each gene in the union
            res = c()
            for (g in 1:length(gen.sig.union)) {
                gene = gen.sig.union[g]
                
                ## Get name, description, biotype                
                gene_name = gen$gene_name[gen$gene_id == gene]
                gene_desc = gen$gene_desc[gen$gene_id == gene]
                gene_type = gen$gene_type[gen$gene_id == gene]
                
                
                ## Get RNAseq and RP log fold changes
                rna_seq             = gen.all$rna_seq[gene]
                ribosomal_profiling = gen.all$ribosomal_profiling[gene]
                
                
                ## Construct regulation string
                rs.reg = constructRegString(gene, rna_seq,                names(gen.sig$rna_seq))
                rp.reg = constructRegString(gene, ribosomal_profiling, names(gen.sig$ribosomal_profiling))
                
                reg = paste0(rs.reg,'/',rp.reg)
                
                res = rbind(res, c(gene, gene_name, gene_desc, gene_type, rna_seq, ribosomal_profiling, reg))
            }
            
            ## Cast result to data.frame
            res = as.data.frame(res, stringsAsFactors=F)
            colnames(res) = c('gene', 'gene_name', 'gene_desc', 'gene_type', 'rna_seq', 'ribosomal_profiling', 'reg')
            res[,c('rna_seq', 'ribosomal_profiling')] = apply(res[,c('rna_seq', 'ribosomal_profiling')], 1:2, as.numeric) 
            
            ## Compute translation efficiency, and sort on it in descending order
            res$te = res$ribosomal_profiling - res$rna_seq
            res    = res[order(res$te, decreasing=T), ]
            
            ## Format gene annotations
            res$gene_name = gsub('_','\\',res$gene_name)
            res$gene_desc = gsub('_','\\',res$gene_desc)
            res$gene_type = gsub('_',' ',res$gene_type)
            
            ## Reorder & rename columns
            res= res[,c('gene', 'gene_name', 'gene_desc', 'gene_type', 'rna_seq', 'ribosomal_profiling', 'te', 'reg')]
			res.csv = res
            colnames(res) = c("Ensembl Gene ID", 
                              "Ensembl Gene Name", 
                              "Ensembl Gene Description", 
                              "Gene Type", 
                              paste0("$Log_{2}$ Fold Change \n (RNA-Seq, ", trtLabl,')'),
                              paste0("$Log_{2}$ Fold Change \n (RP, ", trtLabl,')'),
                              "Translation Efficiency",
                              'Regulation')
			
			## write out to csv
			colnames(res.csv) = c("Ensembl Gene ID", 
					"Ensembl Gene Name", 
					"Ensembl Gene Description", 
					"Gene Type", 
					paste0("Log2 Fold Change (RNA-Seq, ", trtLabl,')'),
					paste0("Log2 Fold Change (RP, ", trtLabl,')'),
					"Translation Efficiency",
					'Regulation')
			out.file = paste0(tbl.dir,'/translation-efficiency-',trtLabl,'-',gsub(' ','',timel),'.csv')
			write.csv(res.csv,file=out.file,row.names=F,quote=T,eol="\n")
			R.utils::gzip(out.file,overwrite=TRUE)
			
                      
            ## Set captions
            caption.short = paste0('Similarities/dissimilarities between sets of identified differentially transcribed and translated genes (',trtLabl,', ',timel,').')
			caption.long  = paste0('Similarities/dissimilarities between sets of identified differentially transcribed and translated genes (',trtLabl,', ',timel,'). Translational efficiency was calculated by taking all genes significantly differentially expressed in both Ribosomal Profiling and RNA-Seq assays, and taking the difference between the Ribosomal Profiling $log_{2}$ fold change and the RNA-Seq $log_{2}$ fold change for each gene. Gene model summaries and annotations are based on Ensembl Version ',ensembl.version,' (August 2017). RP: Ribosomal Profiling, +: significantly up-regulated, -: significantly down-regulated, NS: not significantly regulated, LE: the low-expression filtering cut off was not met.')          
			
			
            ## Add link to ensembl ID
            res[,'Ensembl Gene ID'] = paste0('\\parbox{0.9in}{\\','href{http://useast.ensembl.org/Homo_sapiens/Gene/Summary?g=', res[,'Ensembl Gene ID'],'}{',res[,'Ensembl Gene ID'],'}}')
            
            ## Generate xtable and print
            aln = c('c','p{0.9in}','p{0.6in}','p{2.35in}','p{0.7in}','p{0.7in}','p{0.5in}','p{0.55in}','p{0.55in}')
            lbl = paste0('tab:te',spcFlag,'d',time)
            x = xtable(res, align=aln, type='latex', floating='F', label=lbl, caption=c(caption.long,caption.short))
            
            print(x, tabular.environment='longtable',
                  floating=F,
                  include.rownames=F,
                  hline.after = c(-1,nrow(x)),
                  size='scriptsize',
                  add.to.row = list(pos=list(0),command="\\hline \\endhead "),
                  sanitize.colnames.function=identity,
                  sanitize.text.function=identity,
                  sanitize.rownames.function=identity)
        }
    }
}