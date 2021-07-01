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
# Program:  stage1/Snakefile.sh
# Version:  RPREP 1.0.1
# Author:   Travis L. Jensen, Johannes B. Goll, Sami R. Cherikh, and William F. Hooper
# Purpose:  Run pre-processing for Ribosomal Profiling data or RNA-Seq data
#              Processing steps in brief:
#                1. Directory setup
#                2. Pull reference sequences for masking
#                3. Build Bowtie2 index
#                4. Run fastQC
#                5. Plot raw data read length distribution
#                6. Trim adapters (if rp data)
#                7. FastX quality filter
#                8. Pre-filter read length distribution
#                9. Run length filter [25-35]nt (if rp data)
#                10. Post-filter read length distribution
#                11. Run Bowtie2, extracting unmapped reads
# Input:    {sample}.fastq.gz
# Output:   {sample}_unmapped.fastq.gz
#############################################################################################################

## Import modules
import shutil
import glob
import logging as _logging
import sys

##############################
#        CONFIGURATION         #
##############################

## Specify YAML config file location
configfile: "../preprocess_config.yaml"

## Program locations
CUTADAPT     = config["cutadapt_prog"]
QUALFILTER   = config["fastq_qual_filt_prog"]
BOWTIE2      = config["bowtie_prog"]
BOWTIE2BUILD = BOWTIE2+'-build'
FASTQC       = config["fastqc_prog"]
FASTQDUMP    = config["fastqdump_prog"]
AWS          = config["aws_prog"]
OPENSSL      = config["openssl_prog"]
ENSEMBL_GIT_TOOLS  = config["ensembl_git_tools_prog"]
BIOPERL         = config["bioperl_lib"]

## Assay flag: Skips length filtering if this is rna-seq data
ASSAY        = config["assay"]

## number of cores to use
NCORES  = int(config["ncores"])

## Directories
BASEDIR      = config["pre_dir"]
INPUTDIR     = config["pre_dir"]+'/'+ASSAY+'/input'
SOURCEDIR    = config["srcdir"]
DATADIR        = config["datadir"]

## Use the source dir to import helper modules
sys.path.append(SOURCEDIR+'/python')
import trimadapters
import getfile  
import putfile  
import utils  

## Bowtie2 index name and sequences used 
INDEXNAME    = 'index/masking_idx'
INDEXSEQS    = config["masking_sequences"]

## Adapter sequences
FP_ADAPTERS   = [x.strip() for x in utils.toList(config["fp_adapter_seq"])]
TP_ADAPTERS   = [x.strip() for x in utils.toList(config["tp_adapter_seq"])]

## Ensembl version number
ENSEMBL      = config["ensembl_version"]

## Logging config
LOG_LEVEL = config["log_level"]
LOG_FILE  = config["log_file"]

## Is the data decrypted?
DECRYPT_PASS = config["decrypt_pass"]

## List of samples to process
SAMID_ORIG   = utils.toList(config["samid"]) 
SAMID        = SAMID_ORIG

## List of input files
FASTQ_1      = utils.toList(config["fastq1"])
FASTQ_2      = utils.toList(config["fastq2"])

## run/not run certain steps
RUN_FASTQC = int(config["run_fastqc"])

## Save intermediate files?
REMOVEINTFILES = not config["saveintlocalfiles"]

## Some steps allow paired end samples to be run separately
if (FASTQ_2 != ['']):
    FASTQ_1 += FASTQ_2
    SAMID = [s+'_1' for s in SAMID] + [s+'_2' for s in SAMID]
    ENDS  = ['1','2']
else:
    SAMID = [s+'_1' for s in SAMID]
    ENDS  = ['1']


## Determine whether adapters should be trimmed or not
TRIM_FP = sum([x == 'NA'  for x in FP_ADAPTERS]) == 0
TRIM_TP = sum([x == 'NA'  for x in TP_ADAPTERS]) == 0
TRIM_ADAPTERS_OUTPUT = '.fastq.gz' if (TRIM_FP or TRIM_TP) else '.skipped'


## Configure uploading 
ARCHIVE = config["archive"]
DOARCHIVE = ARCHIVE != 'check_string_for_empty'


## Configure encryption ##TODO: change from test pw
DOENCRYPT    = config["encrypt"]
ENCRYPT_PASS = DECRYPT_PASS

## hash for encryption/decryption & Checksums
HASH = config["hash"]

## Set up logging
_logging.basicConfig(level=LOG_LEVEL, 
                     format='[rprep] %(asctime)s - %(levelname)s - %(message)s', 
                     datefmt='%m/%d/%Y %I:%M:%S %p',
                     handlers=[_logging.FileHandler(LOG_FILE), _logging.StreamHandler(sys.stdout)])
                            


    
####################
# RULE DEFINITIONS #
####################

## On successful workflow completion, merge benchmarks on the rule level
onsuccess:
    workflow_rules = ['fastqc','raw_distribution','trim_adapters','quality_filter',
                    'pre_filt_distribution','length_filter', 'post_filt_distribution',
                    'run_bowtie']
    utils.mergeBenchmarks(samples=SAMID, rules=workflow_rules)
    utils.archiveLog(log=LOG_FILE)
    utils.addResultsToDataDir(datadir=DATADIR, resdir='rlength_distribution')
    
onerror:
    utils.archiveLog(log=LOG_FILE)
    
## Define final output
OUTPUT = [expand('bowtie/{sample}_{pe}_unmapped.fastq.gz', sample=SAMID_ORIG, pe=ENDS),
         expand('rlength_distribution/{sample}_{pe}_raw_counts.RData', sample=SAMID_ORIG, pe=ENDS),
         expand('rlength_distribution/{sample}_{pe}_pre_counts.RData', sample=SAMID_ORIG, pe=ENDS),
         expand('rlength_distribution/{sample}_{pe}_post_counts.RData', sample=SAMID_ORIG, pe=ENDS)]

## Append fastqc done file output if FASTQC
if RUN_FASTQC == 1:
    OUTPUT.append(expand('fastqc/{sample}_{pe}_fastqc.done', sample=SAMID_ORIG, pe=ENDS))

## Define final target (trimmed/filtered fastq files with trnas/rrnas filtered out)
rule all:
    input:
        OUTPUT
    
## Set up directory structure
## Ignores non-zero exit status returned when any directories already exist
rule directory_setup: 
    output: 
        'progress/dirs.done'
    threads:
        1
    priority:
        6
    run:
        shell("mkdir cutadapt rqual_filter rlength_filter rlength_distribution index bowtie fastqc progress benchmark log ../input 2> /dev/null || true")
        utils.touch(output)
        _logging.debug('Created directory structure')



## Set up Ensembl API
rule ensembl_setup:
    output:
        'progress/ensembl_api.done'
    threads:
        1
    priority:
        6
    run:
        ## Try to clone repo (does nothing after running once)
        cmd = 'mkdir -p %s/ensembl; %s --clone api --dir %s/ensembl' % (SOURCEDIR, ENSEMBL_GIT_TOOLS, SOURCEDIR)
        _logging.debug(cmd)
        shell(cmd+utils.returnCode(process='Ensembl API clone', log=LOG_FILE), bench_record=bench_record)
        
        ## Checkout specified ensembl version
        cmd = '%s --checkout --branch release/%s api --dir %s/ensembl' % (ENSEMBL_GIT_TOOLS, ENSEMBL, SOURCEDIR)
        _logging.debug(cmd)
        shell(cmd+utils.returnCode(process='Ensembl API checkout', log=LOG_FILE), bench_record=bench_record)
        utils.touch(output)

    
## Pull reference sequences
rule get_reference_sequences:
    input:
        rules.directory_setup.output,
        rules.ensembl_setup.output
    output:
        fa1='index/ens_trna.fasta',
        fa2='index/ens_xrna.fasta'
    threads: 
        1
    priority:
        6
    benchmark:
        'benchmark/get_reference_sequences.tab'
    run:
        cmd = 'perl %s/perl/get-trna-fastas.pl %s %s %s' % (SOURCEDIR,output.fa1, SOURCEDIR+'/ensembl', BIOPERL)
        _logging.debug(cmd)
        shell(cmd+utils.returnCode(process='tRNA download', log=LOG_FILE), bench_record=bench_record)
        #utils.touch(output.fa1)
        
        cmd = 'Rscript %s/r/get-rrna-fastas.r %s %s' % (SOURCEDIR,output.fa2, ENSEMBL)
        _logging.debug(cmd)
        shell(cmd+utils.returnCode(process='rRNA download', log=LOG_FILE), bench_record=bench_record)
        #utils.touch(output.fa2)
        
    
    
## Build bowtie2 index
rule build_index:
    input:
        rules.get_reference_sequences.output
    output:
        'progress/index_built.done'
    threads: 
        max(1,min(8,NCORES))
    priority:
        1000
    benchmark:
        'benchmark/build_index.tab'
    run:
        index_seqs = ','.join(input + [INDEXSEQS])
        
        cmd = '%s -f %s %s --threads %s' % (BOWTIE2BUILD, index_seqs, INDEXNAME, threads)
        _logging.debug(cmd)
        
        shell(cmd+utils.returnCode(process='Bowtie2 Indexing', log=LOG_FILE),bench_record=bench_record)
        utils.touch(output)
        


## Get input
rule get_file:
    output:
        temp(INPUTDIR+'/{sample}_{pe}.fastq.gz')
    params:
        sample='{sample}_{pe}'
    threads:
        max(1,min(4,NCORES))
    priority:
        1
    run:
        in_file = [FASTQ_1[i[0]] for i in enumerate(SAMID) if SAMID[i[0]] == params.sample]
        getfile.getFile(in_file=in_file[0], out_file=output[0], aws_prog=AWS, fastq_dump_prog=FASTQDUMP, openssl_prog=OPENSSL, pw=DECRYPT_PASS, hash=HASH)
        


## Run FASTQC
rule fastqc: 
    input:
        rules.build_index.output,
        fa=rules.get_file.output
    output:
        temp('fastqc/{sample}_{pe}_fastqc.zip') if REMOVEINTFILES else 'fastqc/{sample}_{pe}_fastqc.zip',
        temp('fastqc/{sample}_{pe}_fastqc.html') if REMOVEINTFILES else 'fastqc/{sample}_{pe}_fastqc.html',
        touch('fastqc/{sample}_{pe}_fastqc.done')
    params:
        sample='{sample}_{pe}'
    threads: 
        1
    priority:
        4
    benchmark:
        'benchmark/{sample}_{pe}_fastqc.tab'
    run:
        cmd = '%s %s -q -o fastqc' % (FASTQC, input.fa)
        _logging.debug(cmd)
        shell(cmd+utils.returnCode(sample='{params.sample}',process='FastQC', log=LOG_FILE),bench_record=bench_record)
        
        if DOENCRYPT: output = utils.encryptFile(file=utils.toList([output[0],output[1]]), openssl=OPENSSL, password=ENCRYPT_PASS, hash=HASH)
        if DOARCHIVE: putfile.upload(file=utils.toList(output), destination=ARCHIVE, cloud='aws', prog=AWS)


        
## Plot raw data read length distribution
rule raw_distribution:
    input:
        fa=INPUTDIR+'/{sample}_{pe}.fastq.gz'
    output:
        rdta='rlength_distribution/{sample}_{pe}_raw_counts.RData',
        pdf='rlength_distribution/{sample}_{pe}_raw_counts.pdf'
    params:
        sample='{sample}_{pe}'
    threads: 
        1
    priority:
        5
    benchmark:
        'benchmark/{sample}_{pe}_raw_distribution.tab'
    run:
        _logging.debug('zcat %s | awk {{if(NR%%4==2) print length($1)}} | Rscript %s/count-read-length-distribution.r %s %s %s raw' % (input.fa, SOURCEDIR, params.sample, output.pdf, output.rdta))
        shell("zcat {input.fa} | awk '{{if(NR%4==2) print length($1)}}' | Rscript {SOURCEDIR}/r/count-read-length-distribution.r {params.sample} {output.pdf} {output.rdta} raw"+utils.returnCode(sample='{params.sample}',process='Raw Read Length Distribution', log=LOG_FILE),bench_record=bench_record)
        
        if DOENCRYPT: output = utils.encryptFile(file=utils.toList(output), openssl=OPENSSL, password=ENCRYPT_PASS, hash=HASH)
        if DOARCHIVE: putfile.upload(file=utils.toList(output), destination=ARCHIVE, cloud='aws', prog=AWS)        
        
        
## Trim adapters with cutadapt
rule trim_adapters:
    input:
        fa=expand('{indir}/{{sample}}_{pe}.fastq.gz',  pe=ENDS, indir=INPUTDIR)
    output:
        temp([x + TRIM_ADAPTERS_OUTPUT for x in expand('cutadapt/{{sample}}_{pe}', pe=ENDS)])
    params:
        sample='{sample}'
    threads: 
        1
    priority:
        5
    benchmark:
        'benchmark/{sample}_trim_adapters.tab'
    run: 
        if (TRIM_FP or TRIM_TP):
            adapters = trimadapters.getAdapters(sam=params.sample, fp_adapters=FP_ADAPTERS, tp_adapters=TP_ADAPTERS, samid_all=SAMID_ORIG)
            trimadapters.trimAdapters(infile=input.fa, outfile=output, adapt5=adapters['fp'], adapt3=adapters['tp'], cutadapt=CUTADAPT)
            
            if DOENCRYPT: output = utils.encryptFile(file=utils.toList(output), openssl=OPENSSL, password=ENCRYPT_PASS, hash=HASH)
            if DOARCHIVE: putfile.upload(file=utils.toList(output), destination=ARCHIVE, cloud='aws', prog=AWS)
        else:
            utils.touch(output)
    
    
    
## Filter out reads whose majority of bases (>50%) have Phred score < 20
rule quality_filter:
    input:
        rules.trim_adapters.output
    output:
        temp('rqual_filter/{sample}_{pe}_qual.fastq.gz')
    params:
        sample='{sample}_{pe}'
    threads: 
        1
    priority:
        5
    benchmark:
        'benchmark/{sample}_{pe}_quality_filter.tab'
    run:
        if (TRIM_FP or TRIM_TP):
            qual_filter_input = input
        else:
            qual_filter_input = INPUTDIR+'/'+params.sample+'.fastq.gz'
        
        cmd = 'zcat %s | %s -Q33 -q 20 -p 50 | gzip -9c > %s' % (qual_filter_input, QUALFILTER, output)
        _logging.debug(cmd)
        shell(cmd+utils.returnCode(sample='{params.sample}', process='FASTQ quality filter', log=LOG_FILE),bench_record=bench_record)
        
        if DOENCRYPT: output = utils.encryptFile(file=utils.toList(output), openssl=OPENSSL, password=ENCRYPT_PASS, hash=HASH)
        if DOARCHIVE: putfile.upload(file=utils.toList(output), destination=ARCHIVE, cloud='aws', prog=AWS)
        
        

## Plot raw data read length distribution
rule pre_filt_distribution:
    input:
        rules.quality_filter.output
    output:
        rdta='rlength_distribution/{sample}_{pe}_pre_counts.RData',
        pdf='rlength_distribution/{sample}_{pe}_pre_counts.pdf'
    params:
        sample='{sample}_{pe}'
    threads: 
        1
    priority:
        5
    benchmark:
        'benchmark/{sample}_{pe}_pre_filt_distribution.tab'
    run:        
        _logging.debug('zcat %s | awk {{if(NR%%4==2) print length($1)}} | Rscript %s/count-read-length-distribution.r %s %s %s pre' % (input, SOURCEDIR, params.sample, output.pdf, output.rdta))
        shell("zcat {input} | awk '{{if(NR%4==2) print length($1)}}' | Rscript {SOURCEDIR}/r/count-read-length-distribution.r {params.sample} {output.pdf} {output.rdta} pre"+utils.returnCode(sample='{params.sample}',process='Pre-filter Read Length Distribution', log=LOG_FILE),bench_record=bench_record)
        
        if DOENCRYPT: output = utils.encryptFile(file=utils.toList(output), openssl=OPENSSL, password=ENCRYPT_PASS, hash=HASH)
        if DOARCHIVE: putfile.upload(file=utils.toList(output), destination=ARCHIVE, cloud='aws', prog=AWS)

        
        
## Filter out reads outside [25,35]nt IF this is RP data
rule length_filter:
    input:
        fa=rules.quality_filter.output
    output:
        temp('rlength_filter/{sample}_{pe}_lenfiltered.fastq.gz')
    params:
        sample='{sample}_{pe}'
    threads: 
        1
    priority:
        5
    benchmark:
        'benchmark/{sample}_{pe}_length_filter.tab'
    run:
        if (ASSAY == 'ribosomal_profiling'):
            _logging.debug('zcat %s | awk BEGIN {{OFS = "\\n"}} {{header = $0; getline seq; getline qheader; getline qseq; if (length(seq) >= 25 && length(seq) <= 35) {{print header, seq, qheader, qseq}}}} | gzip -9c > %s' % (input.fa, output))
            shell("""zcat {input.fa} | awk 'BEGIN {{OFS = "\\n"}} {{header = $0; getline seq; getline qheader; getline qseq; if (length(seq) >= 25 && length(seq) <= 35) {{print header, seq, qheader, qseq}}}}' | gzip -9c > {output}"""+utils.returnCode(sample='{params.sample}',process='Length Filter', log=LOG_FILE),bench_record=bench_record)
            
            if DOENCRYPT: output = utils.encryptFile(file=utils.toList(output), openssl=OPENSSL, password=ENCRYPT_PASS, hash=HASH)
            if DOARCHIVE: putfile.upload(file=utils.toList(output), destination=ARCHIVE, cloud='aws', prog=AWS)
        else:
            _logging.warn('Skipping length filtering step for sample %s' % (params.sample))
            utils.touch(output)
        
        
        
## Plot raw data read length distribution IF this is RP data
rule post_filt_distribution:
    input:
        utils.isPostFilterFile(sample='{sample}_{pe}', assay=ASSAY)
    output:
        rdta='rlength_distribution/{sample}_{pe}_post_counts.RData',
        pdf='rlength_distribution/{sample}_{pe}_post_counts.pdf'
    params:
        sample='{sample}_{pe}'
    threads: 
        1
    priority:
        5
    benchmark:
        'benchmark/{sample}_{pe}_post_filt_distribution.tab'
    run:
        if (ASSAY == 'ribosomal_profiling'):
            _logging.debug('zcat %s | awk {{if(NR%%4==2) print length($1)}} | Rscript %s/count-read-length-distribution.r %s %s %s post' % (input, SOURCEDIR, params.sample, output.pdf, output.rdta))
            shell("zcat {input} | awk '{{if(NR%4==2) print length($1)}}' | Rscript {SOURCEDIR}/r/count-read-length-distribution.r {params.sample} {output.pdf} {output.rdta} post"+utils.returnCode(sample='{params.sample}',process='Post-filter Read Length Distribution', log=LOG_FILE),bench_record=bench_record)
            
            if DOENCRYPT: output = utils.encryptFile(file=utils.toList(output), openssl=OPENSSL, password=ENCRYPT_PASS, hash=HASH)
            if DOARCHIVE: putfile.upload(file=utils.toList(output), destination=ARCHIVE, cloud='aws', prog=AWS)
        else:
            utils.touch(output)

    
    
## Mask rRNA-derived reads with Bowtie2
rule run_bowtie: 
    input:
        utils.isPostFilterFile(sample=expand('{{sample}}_{pe}', pe=ENDS), assay=ASSAY)
    output:
        expand('bowtie/{{sample}}_{pe}_unmapped.fastq.gz', pe=ENDS),
    params:
        sample='{sample}'
    threads:
        max(1,min(8,NCORES))
    priority:
        10
    benchmark:
        'benchmark/{sample}_run_bowtie.tab'
    run:     
        if (len(ENDS) == 1):
            cmd = '%s -p %s -x %s -U %s --un-gz %s -S /dev/null --sensitive-local' % (BOWTIE2, threads, INDEXNAME, input, output)
            _logging.debug(cmd)
        else:
            input_1 = input.fa[0]
            input_2 = input.fa[1]
            cmd = '%s -p %s -x %s -1 %s -2 %s --un-conc-gz bowtie/%s_%%_unmapped.fastq.gz -S /dev/null --sensitive-local' % (BOWTIE2, threads, INDEXNAME, input_1, input_2, params.sample)
            _logging.debug(cmd)

        shell(cmd+utils.returnCode(sample='{params.sample}',process='Bowtie2 alignment', log=LOG_FILE),bench_record=bench_record)
        
        if DOENCRYPT: output = utils.encryptFile(file=utils.toList(output), openssl=OPENSSL, password=ENCRYPT_PASS, hash=HASH)
        if DOARCHIVE: putfile.upload(file=utils.toList(output), destination=ARCHIVE, cloud='aws', prog=AWS)
        
        