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
# Program:  trimadapters.py
# Version:  RPREP 1.1.0
# Author:   Travis L. Jensen, Johannes B. Goll, Sami R. Cherikh, and William F. Hooper
# Purpose:  Cutadapt wrapper script 
# Input:    N/A
# Output:   N/A
#############################################################################################################

import sys
import os
import subprocess
import argparse


## Get the correct adapter sequences for this sample
def getAdapters(sam, fp_adapters, tp_adapters, samid_all):
    fp_adapter = [fp_adapters[i[0]] for i in enumerate(samid_all) if samid_all[i[0]] == sam][0]
    tp_adapter = [tp_adapters[i[0]] for i in enumerate(samid_all) if samid_all[i[0]] == sam][0]
    return({'tp':tp_adapter, 'fp':fp_adapter})


def trimAdapters(cutadapt, infile, outfile, adapt3='NA', adapt5='NA', threads='1'):
    
    ## Determine how many adapters there are, and the ended-ness of the data
    is_paired_end  = False if len(infile) == 1 else True
    is_three_prime = False if adapt3 == 'NA' else True
    is_five_prime  = False if adapt5 == 'NA' else True
    
    ## Check to make sure at least one adapter is present
    if (not is_three_prime and not is_five_prime):
        raise RuntimeError('Missing adapters to use for trimming')
    
    ## If both adapters are present and this is single ended data, use linked trimming
    if (is_three_prime and is_five_prime and not is_paired_end):
        trim_string = '-a '+adapt5+'...'+adapt3+' -o '+outfile[0]+' '+infile[0]
        
    ## If only single ended and 5' adapter, use anchored 5' trimming
    elif (is_five_prime and not is_three_prime and not is_paired_end):
        trim_string = '-g ^'+adapt5+' -o '+outfile[0]+' '+infile[0]
        
    
    ## If only single ended and 3' adapter, use default 3' trimming
    elif (is_three_prime and not is_five_prime and not is_paired_end):
        trim_string = '-a '+adapt3+' -o '+outfile[0]+' '+infile[0]
        
    
    ## If paired end, make sure that both adapters are present
    elif (is_three_prime and is_five_prime and is_paired_end):
        trim_string = '-a '+adapt3 + ' -A '+adapt5+' -o '+outfile[0]+' -p '+outfile[1]+' '+infile[0]+' '+infile[1]
        
    trim_cmd = cutadapt+' '+trim_string+' -m 15 -e 0.1 -O 5 --cores '+str(threads)
    subprocess.run(trim_cmd, shell=True, check=True)