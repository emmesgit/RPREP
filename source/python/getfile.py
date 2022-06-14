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
# Program:  getfile.py
# Version:  RPREP 1.1.0
# Author:   Travis L. Jensen, Johannes B. Goll, Sami R. Cherikh, and William F. Hooper
# Purpose:  get file from local, S3, or SRA
# Input:    N/A
# Output:   N/A
#############################################################################################################

## Helper script to downloaded and/or decrypt files if needed
import sys
import os
import subprocess
import utils


## Download file(s) from amazon s3, return path to result
def downloadS3File(in_file,aws):
    res = []
    for f in in_file:
        cmd = aws+' s3 cp '+f+' '+os.getcwd()+'/'+os.path.basename(f)
        subprocess.run(cmd, shell=True, check=True)
        res.append(os.path.basename(f))
    return(res)



## Construct SRA URL (see https://www.ncbi.nlm.nih.gov/books/NBK158899/)
def getSRAUrl(accession):
    res = 'ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/'+accession[0:6]+'/'+accession+'/'+accession+'.sra'
    return(res)
    
    
    
## Download file(s) from the sequence read archive, return path to result
def downloadSRAFile(accession,fastq_dump):
    res = []
    for f in accession:
    
        ## Download from the NCBI ftp site, extract fastq.gz, remove downloaded file
        ## pefetch is more reliable than fastq-dump to retrieve the sra file
        cmd = 'prefetch -v ' + f + ' -O ' + os.getcwd() + ' && ' + fastq_dump + ' ' + os.getcwd() + '/' + f + '.sra --gzip -O ' + os.getcwd() + ' && rm ' + os.getcwd() + '/' + f + '.sra'
        
        ## If download times out on the first try, give it up to 2 more tries
        for i in range(3):
            try:
                utils.logging_call(cmd, shell=True)
            except subprocess.CalledProcessError as e:
                print('Timeout: ' + e.cmd)
                if (i == 2):
                    raise e
            else:
                break
                
        res.append(os.getcwd() + '/' +f+'.fastq.gz')
    return(res)



def getFile(in_file, out_file, aws_prog, fastq_dump_prog, openssl_prog, pw, hash):
    ## Set flags
    file_downloaded = False
    file_decrypted  = False
    
    ## Convert to list so we can handle multiple files
    in_file = in_file.split(';')
        
    
    ## Download file (if needed)
    if (in_file[0][0:5] == 's3://'):
        in_file = downloadS3File(in_file, aws_prog)
        file_downloaded = True
    elif(in_file[0][0:3] == 'SRR'):
        in_file = downloadSRAFile(in_file, fastq_dump_prog)
        file_downloaded = True

        
    ## Decrypt file (if needed)
    if (in_file[0][-4:] == '.enc'):
        pw = sys.argv[6]
        in_file = utils.decryptFile(in_file, openssl_prog, pw, hash)
        file_decrypted = True
    
    
    ## If we downloaded & decrypted, remove the encrypted file(s)
    if (file_downloaded and file_decrypted):
        [os.remove(s+'.enc') for s in in_file]
        
    
    ## Merge if multiple downloaded files, only keep merged file
    if (len(in_file) > 1 and file_downloaded):
        utils.cat(in_file, out_file)
        [os.remove(s) for s in in_file]
    
    ## Merge if multiple local files, keep unmerged files
    elif (len(in_file) > 1 and not file_downloaded):
        utils.cat(in_file, out_file)
    
    ## Move if single downloaded file
    elif (len(in_file) == 1 and file_downloaded):
        os.rename(in_file[0], out_file)
    
    ## Make a copy if single local file
    elif (len(in_file) == 1 and not file_downloaded):
        cmd = 'cp '+in_file[0]+' '+out_file
        subprocess.run(cmd, shell=True, check=True)
        