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
# Program:  utils.py
# Version:  RPREP 1.0.1
# Author:   Travis L. Jensen, Johannes B. Goll, Sami R. Cherikh, and William F. Hooper
# Purpose:  Python utilities to support prepocessing and snakemake workflow
# Input:    N/A
# Output:   N/A
#############################################################################################################

import datetime
import shutil
import snakemake
import subprocess
import logging
import os
import glob

## Dummy log handler function -- snakemake output gets doubled if updating the logging config 
## without redirecting log output
def logHandler(x):
    return

## Run a process and redirect process output to log
def logging_call(popenargs, **kwargs):
    
    ## Log command
    logging.debug(popenargs)
    
    ## Start process
    process = subprocess.Popen(popenargs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, **kwargs)
        
    ## Check/capture process output    
    def check_io():
        while True:
            output = process.stdout.readline().decode()
            if output:
                logging.info(output)
            else:
                break

    ## keep checking stdout/stderr until the child exits
    while process.poll() is None:
        check_io()
    
    ## Check return code
    if process.returncode != 0:
        raise subprocess.CalledProcessError(returncode=process.returncode, cmd=popenargs)



## Insert object into its own list if it isn't already of class list
def toList(x):
    if isinstance(x, str):
        x = [x]
    elif not isinstance(x, list):
        x = list(x)
    return(x) 
    
    

## cat together files
def cat(in_files, out_file):
    f = ' '.join(in_files)
    cmd = 'cat '+f+' > '+out_file
    subprocess.run(cmd, shell=True, check=True)
    
    
    
## Decrypt file(s), return path to result
def decryptFile(in_file, openssl, password, hash):
    res = []
    for f in in_file:
        cmd = openssl+' aes-256-cbc -md '+hash+' -d -pass pass:'+password+' < '+f+' > ' + os.getcwd() + '/' +os.path.basename(f)[:-4]
        #subprocess.run(cmd, shell=True, check=True)
        logging_call(cmd, shell=True)
        res.append(os.getcwd() + '/' +os.path.basename(f)[:-4])
    return(res)



## Encrypt file(s), return path to result
def encryptFile(file, openssl, password, hash):
    res = []
    
    for f in file:
        cmd = openssl+' aes-256-cbc -md '+hash+' -pass pass:'+password+' < '+f+' > ' + f + '.enc'
        #subprocess.run(cmd, shell=True, check=True)
        logging_call(cmd, shell=True)
        res.append(f + '.enc')
    return(res)
    


## Catch last return code, write to log, pass back to snakemake so 
## that it knows to quit if a non-zero exit status was returned
def returnCode(log, sample='', process=''):
    t = datetime.datetime.now().strftime('%m/%d/%Y %I:%M:%S %p')
    return('; e=$?; echo "[rprep] ' + t + '- INFO - Return Code for ' + sample + ' ' + process + ': $? ">> ' + log + '; exit $e')



## Timestamp and compress the log upon exit
def archiveLog(log):
    t = datetime.datetime.now().strftime('%Y-%m-%d.h%H-m%M-s%S')
    new_log = log+'.'+t
    try:
        shutil.copy2(log, new_log)
        snakemake.shell('gzip -9 '+new_log)
    except (FileNotFoundError, PermissionError) as e:
        pass



## If this is RNA-Seq data, we're not filtering by length
## This function is used as input to bowtie/post-filter plot script
## Essentially skips filtering-related nodes in the DAG 
def isPostFilterFile(sample, assay):
    if (assay == 'ribosomal_profiling'):
        return(['rlength_filter/'+s+'_lenfiltered.fastq.gz' for s in toList(sample)])
    elif (assay == 'rna_seq'):
        return(['rqual_filter/'+s+'_qual.fastq.gz' for s in toList(sample)])



## Shorthand for touch - takes list<str> or str
def touch(file):
    if type(file) is list:
        [snakemake.shell('touch ' + x) for x in file]
    else:
        snakemake.shell('touch {file}')
    
    
    
## Merge benchmarking output by rule at the end of a succesful run
def mergeBenchmarks(samples, rules):

    ## Keep merged and unmerged benchmarks in seperate dirs
    snakemake.shell('mkdir benchmark/merged benchmark/unmerged -p 2> /dev/null')
    
    ## If there's a single sample, convert to list to avoid looping by char
    if (isinstance(samples, str)):
        samples = [samples]
    
    for rule in rules:
        f = 'benchmark/merged/'+rule+'.tab'
        
        ## Don't add the header if we don't have to
        if (not os.path.isfile(f)):
            snakemake.shell('echo -e "s\th:m:s\tmax_rss\tmax_vms\tmax_uss\tmax_pss\tio_in\tio_out\tmean_load" >> {f}')
            
        for sample in samples:
            snakemake.shell('tail -n1 benchmark/{sample}_*_{rule}.tab >> {f}  2> /dev/null || true')
            snakemake.shell('mv benchmark/{sample}_*_{rule}.tab benchmark/unmerged 2> /dev/null|| true ')
            
    ## Move non-sample dependent rule benchmarks to merged
    snakemake.shell('mv benchmark/*tab benchmark/merged 2> /dev/null || true')



## Copy length counts to the data directory
def addResultsToDataDir(datadir, resdir):
    destdir = datadir + '/' + resdir 
    if (not os.path.isdir(destdir)):
        os.mkdir(destdir)
    [shutil.copy2(x, destdir) for x in glob.glob(resdir+'/*_counts.RData')]
    
    
    
## Determine the outputs from HISAT2 -- changes depending on whether we're remapping or not
def hisatOutput(ends, iontorrent):
    res = {'aln' : snakemake.io.temp('tmp/{sample}.sam')}
            
    if (iontorrent):
        res['un_aln'] = snakemake.io.temp(expand('tmp/{{sample}}_iontorrent_{pe}.fastq.gz', pe=ends))
        
    return(res)
    


## Determine bam output based on whether we're remapping/compressing
def sortBamOutput(iontorrent, cram):
    if (iontorrent and cram):
        return(snakemake.io.temp('tmp/{sample}_sorted.bam'))
    else:
        if (cram == 1):
            return(snakemake.io.temp('bam/{sample}.bam'))
        else: 
            return('bam/{sample}.bam')



## Merge RSeQC results, return paths to merged results
def mergeRSEQC(srcdir, run_read_dist):
    retvals = ['rseqc/bam_qc_parsed.tab', 'rseqc/bam_gc_parsed.tab', 'rseqc/bam_jc_parsed.tab', 'feature_counts/fragment_count_matrix.tab.gz','feature_counts/gene_lengths.tab']
    
    ## QC
    snakemake.shell("find rseqc | grep -P 'bam_qc.txt$' > bam_qc_outfiles.txt")
    snakemake.shell("perl {srcdir}/perl/parse-rseqc-bam-qc-results.pl bam_qc_outfiles.txt > rseqc/bam_qc_parsed.tab")
    snakemake.shell("rm bam_qc_outfiles.txt")

    ## GC
    snakemake.shell("find rseqc | grep -P 'bam_gc.txt$' > bam_gc_outfiles.txt")
    snakemake.shell("perl {srcdir}/perl/parse-rseqc-bam-gc-results.pl bam_gc_outfiles.txt > rseqc/bam_gc_parsed.tab")
    snakemake.shell("rm bam_gc_outfiles.txt")

    ## JC
    snakemake.shell("find rseqc | grep -P 'bam_jc.txt$' > bam_jc_outfiles.txt")
    snakemake.shell("perl {srcdir}/perl/parse-rseqc-bam-jc-results.pl bam_jc_outfiles.txt > rseqc/bam_jc_parsed.tab")
    snakemake.shell("rm bam_jc_outfiles.txt")
    
    ## RC
    if (run_read_dist) == 1:
        snakemake.shell("find rseqc | grep -P 'bam_rc.txt$' > bam_rc_outfiles.txt")
        snakemake.shell("perl {srcdir}/perl/parse-rseqc-read-distribution-results.pl bam_rc_outfiles.txt > rseqc/bam_rc_parsed.tab")
        snakemake.shell("rm bam_rc_outfiles.txt")
        retvals.append('rseqc/bam_rc_parsed.tab')

    ## featureCounts
    snakemake.shell("find feature_counts | grep -P '_count.tab$' > feature_counts_outfiles.txt")
    snakemake.shell("perl {srcdir}/perl/parse-read-count-matrix-subread-1.4.6.pl feature_counts_outfiles.txt > feature_counts/fragment_count_matrix.tab")
    snakemake.shell("gzip -f feature_counts/fragment_count_matrix.tab")
    
    ## Export gene lengths from a single featurecounts file 
    snakemake.shell("head -1 feature_counts_outfiles.txt | xargs tail -n +2 | awk '{{print $1,$6}}' > feature_counts/gene_lengths.tab")
    snakemake.shell("rm feature_counts_outfiles.txt")
    
    return(retvals)
