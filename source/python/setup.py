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
# Program:  setup.py
# Version:  RPREP 1.1.0
# Author:   Travis L. Jensen, Johannes B. Goll, Sami R. Cherikh, and William F. Hooper
# Purpose:  Pulls latest analysis components from the RSEQREP Github repo
# Input:    N/A
# Output:   N/A
#############################################################################################################

import sys
import os
import shutil
import subprocess
import json
import logging

GIT_URL = 'https://api.github.com/repos/emmesgit/RSEQREP/releases/latest'


## Pull latest RSEQREP source code
def getLatestRSEQREPSource(dir, git_url=GIT_URL):
    
    cmd = 'mkdir '+dir+'; curl -s ' + git_url + ' | grep "tarball_url" | cut -d : -f 2,3 | tr -d \\" | sed "s/,//g" | wget -q -i - -O '+dir+'/rseqrep --no-check-certificate'
    
    try:
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        logging.debug(e.cmd)
        logging.error("Couldn't download RSEQREP analytical components")
        exit(1)
        
        

## Pull latest RSEQREP metadata
def getLatestRSEQREPMetadata(dir, git_url=GIT_URL):
    git_url = 'https://api.github.com/repos/emmesgit/RSEQREP/releases/latest'
    json_file = dir+'/rseqrep-release-metadata.json'
    version_file = dir+'/.version'
    cmd = 'curl -s ' + git_url + '> '+ json_file
    
    ## Try to download metadata
    try:
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError:
        logging.error("Couldn't download RSEQREP metadata")
        exit(1)
        
    ## Pull version from metadata
    version = json.load(open(json_file))['tag_name']
    
    with open(version_file, 'w') as f:
        f.write(version)
    



def setup(rprep_dir):
            
    ## Get working directories
    rprep_r_dir   = rprep_dir + '/source/r'
    rseqrep_dir   = rprep_dir+'/source/rseqrep'
    rseqrep_r_dir = rseqrep_dir+'/source/r/'
    
    logging.info('Attemping to download and unpack RSEQREP analytical components')
    
    ## Pull the latest RSEQREP release and its metadata
    getLatestRSEQREPSource(rseqrep_dir)
    getLatestRSEQREPMetadata(rseqrep_dir)
    
    ## Untar RSEQREP
    untar_dir = rseqrep_dir
    cmd = 'tar -xzf '+rseqrep_dir+'/rseqrep --strip-components 1 -C '+ untar_dir
    try:
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        logging.debug(e.cmd)
        logging.error("Couldn't unpack RSEQREP analytical components")
        exit(1)
    
    ## Remove RSEQREP tarball
    os.remove(rseqrep_dir+'/rseqrep')
    
    logging.info('RSEQREP analytical components successfully downloaded and unpacked')
    logging.info('Attempting to patch RSEQREP analytical components')
    
    
    ## List directories of files to patch and their replacements
    old_dirs  = [rseqrep_r_dir, rseqrep_dir, rseqrep_r_dir+'/04-sdeg-organization-known-modules', rseqrep_r_dir+'02-sdeg-identification', rseqrep_r_dir+'03-sdeg-clusters']
    new_files = [rprep_r_dir + '/init-analysis.r', rprep_dir + '/source/shell/run-analysis.sh', rprep_r_dir+'/init-gsea-sampling.r',rprep_r_dir+'/init-edgeR-glm-model.r',rprep_r_dir+'/init-pvclusters.r']
    
    try:
        [shutil.copy2(new, old) for (new, old) in zip(new_files, old_dirs)]
    except FileNotFoundError:
        logging.error("Couldn't patch RSEQREP analytical components")
        exit(1)
    
    logging.info('RSEQREP analytical components successfully patched')
    
    return(True)