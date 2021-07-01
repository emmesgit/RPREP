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
# Program:  stage2/Snakefile.sh
# Version:  RPREP 1.0.1
# Author:   Travis L. Jensen, Johannes B. Goll, Sami R. Cherikh, and William F. Hooper
# Purpose:  Map reads to reference genome, gather QC metrics, run featureCounts
# Input:    {sample}.<fastq.gz>
# Output:   {sample}.bam
#           rseqc/bam_<qc|jc|rc|gc>_parsed.tab
#############################################################################################################

## Import modules
import shutil
import logging as _logging

##############################
#        CONFIGURATION         #
##############################

## Specify YAML config file location
configfile: "../preprocess_config.yaml"

## Assay type
ASSAY          = config["assay"]

## number of cores to use
NCORES  = int(config["ncores"])

## Program locations
RSEQC         = config["rseqc_dir"]
SAMTOOLS      = config["samtools_prog"]
HISAT2        = config["hisat_prog"]
HISAT2BUILD   = HISAT2 + '-build'
HISAT2SPLICE  = HISAT2 + '_extract_splice_sites.py'
BOWTIE2       = config["bowtie_prog"]
BOWTIE2BUILD  = BOWTIE2+'-build'
FEATURECOUNTS = config["fcts_prog"]

## Directories
INPUTDIR      = config["stage1dir"]
SOURCEDIR     = config["srcdir"]
DATADIR       = config["datadir"]
RESDIR        = 'bias_plots'

## Use the source dir to import helper modules
sys.path.append(SOURCEDIR+'/python')
import getfile  
import putfile  
import utils  

## Ensembl version number
ENSEMBL = config["ensembl_version"]

## Logging config
LOG_LEVEL = config["log_level"]
LOG_FILE  = config["log_file"]

## run/not run certain steps
RUN_READ_DIST = int(config["run_read_dist"])

## Save intermediate files?
REMOVEINTFILES = not config["saveintlocalfiles"]

## Reference datasets
INDEXSEQ = 'genome/Homo_sapiens.ensembl.version'+str(ENSEMBL)+'.genome.fa'
ANNOTATIONS_BED = 'annot/Homo_sapiens.ensembl.version'+str(ENSEMBL)+'.chr.bed'
ANNOTATIONS_GTF = 'annot/Homo_sapiens.ensembl.version'+str(ENSEMBL)+'.chr.gtf'

## Remap reads? 
IONTORRENT   = config["iontorrent"]

## Generate CRAM files?
CRAM         = int(config["save_cram"]) 
CRAMFLAG     = ['bam','cram'][CRAM]

## Perform gene-body-coverage?
GENEBODYCOVERAGE = config["genebodycoverage"]

## Save intermediate files?
REMOVEINTFILES = not config["saveintlocalfiles"]

## List of samples to process
SAMID        = utils.toList(config["samid"]) 

## List of input files
FASTQ_1      = utils.toList(config["fastq1"])
FASTQ_2      = utils.toList(config["fastq2"])

## Check to see if reads are paired-end
ENDS = ['1','2'] if (FASTQ_2 != ['']) else ['1']

## Configure uploading 
ARCHIVE = config["archive"]
DOARCHIVE = ARCHIVE != 'check_string_for_empty'


## Configure encryption
DECRYPT_PASS = config["decrypt_pass"]
DOENCRYPT    = config["encrypt"]
ENCRYPT_PASS = DECRYPT_PASS

## hash for encryption/decryption & Checksums
HASH = config["hash"]

## Function to copy gene_bias_plots output to results data dir
def addResultsToDataDir(datadir=DATADIR, resdir=RESDIR):
    destdir = datadir + '/' + resdir 
    if (not os.path.isdir(destdir)):
        os.mkdir(destdir)
    [shutil.copy2(x, destdir) for x in glob.glob(resdir+'/*geneBodyCoverage.txt')]


## Set up logging
_logging.basicConfig(level=LOG_LEVEL, 
                    format='[rprep] %(asctime)s - %(levelname)s - %(message)s', 
                    datefmt='%m/%d/%Y %I:%M:%S %p',
                    handlers=[_logging.FileHandler(LOG_FILE),
                             _logging.StreamHandler(sys.stdout)])


## Define final output
OUTPUT = [expand('rseqc/{sample}_bam_qc.txt', sample=SAMID),
        expand('rseqc/{sample}_bam_gc.txt', sample=SAMID),
        expand('rseqc/{sample}_bam_jc.txt', sample=SAMID),
        expand('feature_counts/{sample}_count.tab', sample=SAMID)]

## Append gene_bias_plots output if GENEBODYCOVERAGE
if GENEBODYCOVERAGE:
    OUTPUT.append(expand(RESDIR+'/{sample}.geneBodyCoverage.txt', sample=SAMID))
if CRAM == 1 and REMOVEINTFILES:
    OUTPUT.append(expand('progress/{sample}_cram.done', sample=SAMID))
if not REMOVEINTFILES:
    OUTPUT.append(expand(CRAMFLAG+'/{sample}.'+CRAMFLAG, sample=SAMID))
if RUN_READ_DIST == 1:
    OUTPUT.append(expand('rseqc/{sample}_bam_rc.txt', sample=SAMID))


####################
# RULE DEFINITIONS #
####################

## On successful completion of the workflow, merge rseqc results and benchmarks, move to data dir
onsuccess:
    workflow_rules = ['build_hisat_index','build_bowtie_index','run_hisat','run_bowtie',
                    'sam_to_bam','sort_bam','merge_bam','bam_qc','bam_gc','bam_jc','read_distribution',
                    'feature_counts','index_bam','gene_bias_plots']
   
    ## merge featurecounts, rseqc results
    merged_results = utils.mergeRSEQC(SOURCEDIR, RUN_READ_DIST)
    [shutil.copy2(x, DATADIR) for x in merged_results]
    
    if GENEBODYCOVERAGE:
        addResultsToDataDir()
    
    utils.mergeBenchmarks(samples=SAMID, rules=workflow_rules)
    utils.archiveLog(log=LOG_FILE)
onerror:
    utils.archiveLog(log=LOG_FILE)
  
## Define final output
rule all:
    input:
        OUTPUT
    
## Set up directory structure
## Ignores non-zero exit status returned when any directories already exist
rule directory_setup: 
    output: 
        'progress/dirs.done'
    threads: 1
    priority: 2
    run:
        cmd = "mkdir cram annot genome index input bam feature_counts rseqc tmp -p 2> /dev/null"
        _logging.debug(cmd)
        shell(cmd+utils.returnCode(process='Directory setup', log=LOG_FILE), bench_record=bench_record)
        utils.touch(output)

    
    
## Build HISAT2 index
rule build_hisat_index:
    input:
        rules.directory_setup.output
    output:
        'progress/hisat2_index_built.done'
    benchmark:
        'benchmark/build_index.tab'
    threads: max(1,min(8,NCORES))
    priority: 1000
    run:
        cmd = '%s -p %s %s index/hisat2_genome_index && %s %s > annot/splicesites.txt' % (HISAT2BUILD, threads, INDEXSEQ, HISAT2SPLICE, ANNOTATIONS_GTF)
        _logging.debug(cmd)
        shell(cmd+utils.returnCode(process='HISAT2 Build', log=LOG_FILE), bench_record=bench_record)    
        utils.touch(output)



## Build Bowtie2 index
rule build_bowtie_index:
    input:
        rules.directory_setup.output
    output:
        'progress/bowtie2_index_built.done'
    benchmark:
        'benchmark/bowtie2_index.tab'
    threads: max(1,min(8,NCORES))
    priority: 33
    run:
        cmd = '%s -f %s index/bowtie2_genome_index --threads %s' % (BOWTIE2BUILD, INDEXSEQ, threads))
        _logging.debug(cmd)
        shell(cmd+utils.returnCode(process='Bowtie2 Build', log=LOG_FILE), bench_record=bench_record)
        utils.touch(output)



## Index reference genome with faidx (only necessary if CRAM compression is enabled)
rule faidx: 
    output:
        'progress/faidx.done'
    benchmark:
        'benchmark/faidx.tab'
    threads: 1
    priority: 4
    run:
        cmd = '%s faidx %s' % (SAMTOOLS, INDEXSEQ)
        _logging.debug(cmd)
        shell(cmd+utils.returnCode(process='Bowtie2 Build', log=LOG_FILE), bench_record=bench_record)
        utils.touch(output)
        

## Map reads to the reference genome using HISAT2 
rule run_hisat:
    input:
        tch=rules.build_hisat_index.output,
        fa=expand('{input}/bowtie/{{sample}}_{pe}_unmapped.fastq.gz', pe=ENDS, input=INPUTDIR)
    output:
        **utils.hisatOutput(ends=ENDS, iontorrent=IONTORRENT)
    benchmark:
        'benchmark/{sample}_run_hisat.tab'
    threads: max(1,min(8,NCORES-4))
    priority: 5
    params:
        sample='{sample}'
    run:
        ## Construct input based on read ended-ness
        if (len(ENDS) == 1):
            in_fa_str = "-U "+input.fa[0]
            
            if (IONTORRENT):
                unmapped_str = "--un-gz tmp/"+params.sample+"_iontorrent_1.fastq.gz"
            else:
                unmapped_str = ''
            
        elif (len(ENDS) == 2):
            in_fa_str = '-1 '+input.fa[0]+' -2 '+input.fa[1]
            
            if (IONTORRENT):
                unmapped_str = "--un-conc-gz tmp/"+params.sample+"_iontorrent_%.fastq.gz"
            else:
                unmapped_str = ''
        
        cmd = '%s -p %s -x index/hisat2_genome_index %s %s --known-splicesite-infile annot/splicesites.txt -S %s' % (HISAT2, threads, in_fa_str, unmapped_str, output.aln)
        _logging.debug(cmd)
        shell(cmd+utils.returnCode(process='HISAT2 Align', sample='{params.sample}', log=LOG_FILE), bench_record=bench_record)
        
        if DOENCRYPT: output = utils.encryptFile(file=utils.toList(output), openssl=OPENSSL, password=ENCRYPT_PASS, hash=HASH)
        if DOARCHIVE: putfile.upload(file=utils.toList(output), destination=ARCHIVE, cloud='aws', prog=AWS)        


rule run_bowtie:
    input:
        rules.build_bowtie_index.output,
        fa=expand('tmp/{{sample}}_iontorrent_{pe}.fastq.gz', pe=ENDS)
    output:
        temp('tmp/{sample}_iontorrent.bam')
    benchmark:
        'benchmark/{sample}_run_bowtie.tab'
    threads: max(1,min(8,NCORES-4))
    priority: 5
    params:
        sample='{sample}'
    run:
        ## Construct input based on read ended-ness
        if (len(ENDS) == 1):
            in_fa_str = "-U "+input.fa[0]
        elif (len(ENDS) == 2):
            in_fa_str = '-1 '+input.fa[0]+' -2 '+input.fa[1]
        
        cmd = '%s --threads %s -x index/bowtie2_genome_index %s --sensitive-local | %s view -h -u -b | %s sort -@ %s - > %s' % (BOWTIE2, threads, in_fa_str, SAMTOOLS, SAMTOOLS, threads, output)
        _logging.debug(cmd)
        shell(cmd+utils.returnCode(process='Bowtie2 Remapping', sample='{params.sample}', log=LOG_FILE), bench_record=bench_record)
        
        if DOENCRYPT: output = utils.encryptFile(file=utils.toList(output), openssl=OPENSSL, password=ENCRYPT_PASS, hash=HASH)
        if DOARCHIVE: putfile.upload(file=utils.toList(output), destination=ARCHIVE, cloud='aws', prog=AWS)        


## Convert SAM to BAM
rule sam_to_bam:
    input:
        rules.run_hisat.output.aln
    output:
        temp('tmp/{sample}.bam')
    benchmark:
        'benchmark/{sample}_sam_to_bam.tab'
    threads: max(1,min(8,NCORES-4))
    priority: 5
    params:
        sample='{sample}'
    run:
        cmd = '%s view -@ %s -Sbh %s > %s' % (SAMTOOLS, threads, input, output)
        _logging.debug(cmd)
        shell(cmd+utils.returnCode(process='SAM to BAM', sample='{params.sample}', log=LOG_FILE), bench_record=bench_record)
        
        if DOENCRYPT: output = utils.encryptFile(file=utils.toList(output), openssl=OPENSSL, password=ENCRYPT_PASS, hash=HASH)
        if DOARCHIVE: putfile.upload(file=utils.toList(output), destination=ARCHIVE, cloud='aws', prog=AWS)
        


## Sort BAM
rule sort_bam:
    input:
        bam=rules.sam_to_bam.output,
    output:
        utils.sortBamOutput(iontorrent=IONTORRENT, cram=CRAM)
    benchmark:
        'benchmark/{sample}_sort_bam.tab'
    threads: max(1,min(8,NCORES-4))
    priority: 5
    params:
        sample='{sample}'
    run:
        cmd = '%s sort -@ %s %s > %s' % (SAMTOOLS, threads, input.bam, output)
        _logging.debug(cmd)
        shell(cmd+utils.returnCode(process='Sort BAM', sample='{params.sample}', log=LOG_FILE), bench_record=bench_record)
        
        if DOENCRYPT: output = utils.encryptFile(file=utils.toList(output), openssl=OPENSSL, password=ENCRYPT_PASS, hash=HASH)
        if DOARCHIVE: putfile.upload(file=utils.toList(output), destination=ARCHIVE, cloud='aws', prog=AWS)        



## Merge sorted alignments if we re-mapped with bowtie
rule merge_bam:
    input:
        bowtie=rules.run_bowtie.output,
        hisat=rules.sort_bam.output
    output:
        temp('bam/{sample}.bam') if REMOVEINTFILES else  'bam/{sample}.bam'
    benchmark:
        'benchmark/{sample}_merge_bam.tab'
    threads: 1
    priority: 5
    params:
        sample='{sample}'
    run:
        cmd = '%s merge %s %s %s' % (SAMTOOLS, output, input.hisat, input.bowtie)
        _logging.debug(cmd)
        shell(cmd+utils.returnCode(process='Merge BAM', sample='{params.sample}', log=LOG_FILE), bench_record=bench_record)
        
        if DOENCRYPT: output = utils.encryptFile(file=utils.toList(output), openssl=OPENSSL, password=ENCRYPT_PASS, hash=HASH)
        if DOARCHIVE: putfile.upload(file=utils.toList(output), destination=ARCHIVE, cloud='aws', prog=AWS)
        
        

rule bam_to_cram:
    input:
        bam='bam/{sample}.bam',
        faidx=rules.faidx.output
    output:
        cram=temp('cram/{sample}.cram') if REMOVEINTFILES else 'cram/{sample}.cram',
        prog=touch('progress/{sample}_cram.done')
    benchmark: 
        'benchmark/{sample}_bam_to_cram.tab'
    threads: max(1,min(8,NCORES-4))
    priority: 5
    params:
        sample='{sample}'
    run:
        cmd = '%s view -C -@ %s -T %s -o %s %s' % (SAMTOOLS, threads, INDEXSEQ, output.cram, input.bam)
        _logging.debug(cmd)
        shell(cmd+utils.returnCode(process='CRAM Compression', sample='{params.sample}', log=LOG_FILE), bench_record=bench_record)
        
        if DOENCRYPT: output = utils.encryptFile(file=utils.toList(output), openssl=OPENSSL, password=ENCRYPT_PASS, hash=HASH)
        if DOARCHIVE: putfile.upload(file=utils.toList(output), destination=ARCHIVE, cloud='aws', prog=AWS)
        if DOARCHIVE: putfile.upload(file=utils.toList(INDEXSEQ), destination=ARCHIVE, cloud='aws', prog=AWS)


## Run RSEQC bam_stat.py
rule bam_qc:
    input:
        'bam/{sample}.bam'
    output:
        'rseqc/{sample}_bam_qc.txt'
    benchmark:
        'benchmark/{sample}_bam_qc.tab'
    threads: 1
    priority: 5
    params:
        sample='{sample}'
    run:
        cmd = '%s/bam_stat.py -i %s > %s' % (RSEQC, input, output)
        _logging.debug(cmd)
        shell(cmd+utils.returnCode(process='BAM QC', sample='{params.sample}', log=LOG_FILE), bench_record=bench_record)    
        
        if DOENCRYPT: output = utils.encryptFile(file=utils.toList(output), openssl=OPENSSL, password=ENCRYPT_PASS, hash=HASH)
        if DOARCHIVE: putfile.upload(file=utils.toList(output), destination=ARCHIVE, cloud='aws', prog=AWS)        
    
    
            
## Run RSEQC read_gc.py 
rule bam_gc:
    input:
        'bam/{sample}.bam'
    output:
        r='rseqc/{sample}.GC_plot.r' if REMOVEINTFILES else 'rseqc/{sample}.GC_plot.r',
        txt='rseqc/{sample}_bam_gc.txt'
    benchmark:
        'benchmark/{sample}_bam_gc.tab'
    threads: 1
    priority: 5
    params:
        sample='{sample}'
    run:
        ## Run RSEQC read_GC
        cmd = '%s/read_GC.py -i %s -o rseqc/%s' % (RSEQC, input, params.sample)
        _logging.debug(cmd)
        shell(cmd+utils.returnCode(process='BAM GC', sample='{params.sample}', log=LOG_FILE), bench_record=bench_record)    
        
        ## Add some extra code to output R script, run it
        cmd = """echo "out=as.vector(summary(gc));dta = data.frame('%s',out[1],out[2],out[3],out[4],out[5],out[6]);write.table(dta,file='%s',sep='\t',row.names=F,col.names=F,quote=F);" >> %s""" % (params.sample, output.txt, output.r)
        _logging.debug(cmd)
        shell(cmd)
        shell('Rscript --vanilla --quiet {output.r}'+utils.returnCode(process='BAM GC Rscript', sample='{params.sample}', log=LOG_FILE))
        
        if DOENCRYPT: output = utils.encryptFile(file=utils.toList(output), openssl=OPENSSL, password=ENCRYPT_PASS, hash=HASH)
        if DOARCHIVE: putfile.upload(file=utils.toList(output), destination=ARCHIVE, cloud='aws', prog=AWS)
        
        ## Remove intermediate read_GC files after s3 upload
        if REMOVEINTFILES:
            cmd = 'rm rseqc/%s.GC.xls rseqc/%s.GC_plot.pdf' % (params.sample, params.sample)
            utils.logging_call(cmd, shell=True)

    
## Run RSEQC junction_annotation.py 
rule bam_jc:
    input:
        'bam/{sample}.bam'
    output:
        'rseqc/{sample}_bam_jc.txt'
    benchmark:
        'benchmark/{sample}_bam_jc.tab'
    threads: 1
    priority: 5
    params:
        sample='{sample}',
    run:
        cmd = '%s/junction_annotation.py -i %s -o rseqc/%s -r %s 2> %s' % (RSEQC, input, params.sample, ANNOTATIONS_BED, output)
        _logging.debug(cmd)
        shell(cmd+utils.returnCode(process='BAM JC', sample='{params.sample}', log=LOG_FILE),bench_record=bench_record)    
        
        if DOENCRYPT: output = utils.encryptFile(file=utils.toList(output), openssl=OPENSSL, password=ENCRYPT_PASS, hash=HASH)
        if DOARCHIVE: putfile.upload(file=utils.toList(output), destination=ARCHIVE, cloud='aws', prog=AWS)
        
        

## Run RSEQC junction_annotation.py 
rule read_distribution:
    input:
        'bam/{sample}.bam'
    output:
        'rseqc/{sample}_bam_rc.txt'
    benchmark:
        'benchmark/{sample}_bam_rc.tab'
    threads: 1
    priority: 5
    params:
        sample='{sample}',
    run:
        cmd = '%s/read_distribution.py -i %s -r %s > %s' % (RSEQC, input, ANNOTATIONS_BED, output)
        _logging.debug(cmd)
        shell(cmd+utils.returnCode(process='BAM Read Distribution', sample='{params.sample}', log=LOG_FILE),bench_record=bench_record)    
        
        if DOENCRYPT: output = utils.encryptFile(file=utils.toList(output), openssl=OPENSSL, password=ENCRYPT_PASS, hash=HASH)
        if DOARCHIVE: putfile.upload(file=utils.toList(output), destination=ARCHIVE, cloud='aws', prog=AWS)
        
        
        
## Run featureCounts
rule feature_counts:
    input:
        'bam/{sample}.bam'
    output:
        'feature_counts/{sample}_count.tab'
    benchmark:
        'benchmark/{sample}_feature_counts.tab'
    threads: max(1,min(8,NCORES-4))
    priority: 5
    params:
        sample='{sample}',
    run:
        cmd = '%s -T %s -a %s -o %s %s' % (FEATURECOUNTS, threads, ANNOTATIONS_GTF, output, input)
        _logging.debug(cmd)
        shell(cmd+utils.returnCode(process='FeatureCounts', sample='{params.sample}', log=LOG_FILE),bench_record=bench_record)
        
        if DOENCRYPT: output = utils.encryptFile(file=utils.toList(output), openssl=OPENSSL, password=ENCRYPT_PASS, hash=HASH)
        if DOARCHIVE: putfile.upload(file=utils.toList(output), destination=ARCHIVE, cloud='aws', prog=AWS)
        


## Index bam
rule index_bam:
    input:
        'bam/{sample}.bam'
    output:
        temp('bam/{sample}.bam.bai')
    threads: max(1,min(8,NCORES-4))
    priority: 6
    params:
        sample='{sample}'
    run:
        ## indexing bam
        cmd = '%s index -@ %s %s' % (SAMTOOLS, threads, input)
        _logging.debug(cmd)
        shell(cmd+utils.returnCode(process='Samtools index', sample='{params.sample}', log=LOG_FILE), bench_record=bench_record)



## Run RSeQC geneBodyCoverage
rule gene_bias_plots:
    input:
        bam='bam/{sample}.bam',
        bai=rules.index_bam.output
    output:
        RESDIR+'/{sample}.geneBodyCoverage.txt'
    threads: 1
    priority: 7
    params:
        sample='{sample}'
    run:
        cmd = '%s/geneBody_coverage.py -r %s -i %s -o %s/%s' % (RSEQC, ANNOTATIONS_BED, input.bam, RESDIR, params.sample)
        _logging.debug(cmd)
        shell(cmd+utils.returnCode(process='RSeQC geneBodyCoverage', sample='{params.sample}', log=LOG_FILE), bench_record=bench_record)

        if DOENCRYPT: output = utils.encryptFile(file=utils.toList(output), openssl=OPENSSL, password=ENCRYPT_PASS, hash=HASH)
        if DOARCHIVE: putfile.upload(file=utils.toList(output), destination=ARCHIVE, cloud='aws', prog=AWS)
        

    