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
# Program:  parse-rseqc-bam-qc-results.pl 
# Version:  RPREP 1.0.1
# Author:   Travis L. Jensen, Johannes B. Goll, Sami R. Cherikh, and William F. Hooper
# Purpose:  Parse RSeqC bam_stat.py results
# Input:    list of absolute file paths
# Output:   N/A
#############################################################################################################

use strict;
use warnings;
use File::Basename;

my $isFirst=1;

my $header = "sample_id\ttotal\tqc_failed\tpcr_duplicates\tnon_primary\tunmapped\tnon_unique\tunique\tread_1\tread_2\tplus_strand\tminus_strand\tnon_spliced\tspliced\tproper_pairs\tnon_proper_pairs";

while(<>) {
    chomp;
    my $file = $_;
    my $sampleId = basename($file);
    $sampleId =~ s/_bam_qc\.txt//;
    
    open(DATA, "<$file") or die "Couldn't open file $file, $!";
    my $resLine=undef;
    my $hedLine=undef;
    
    while(<DATA>){
        chomp;
        
        my $line = $_;
        $line =~ s/Non primary hits/Non primary hits:/;
        
        if($line=~ m/Segmentation fault/) {
            print "$sampleId\tSegmentation fault (core dumped)\n";
            last;
        }
        
        unless($line =~ m/^$/ || $line=~ m/Load/ || $line=~ m/^\s/ || $line=~ m/^#/ ) {
            
            my @linePart = split('\:',$line,2);
            
            my $value= trim($linePart[1]);
            
            if(!defined $value) { $value =0;}; 
            
            if(! defined $resLine) {
                $resLine="$sampleId\t$value";
            } else{
                $resLine="$resLine\t$value";
            }
        }
    }
    if($isFirst) {
        print $header."\n";
        $isFirst=0;
    } 
    
    if(defined $resLine) {
        print $resLine."\n";
    }
    close DATA;
    
}

sub trim {
    my $s = shift;
    $s =~ s/^\s+|\s+$//g ; 
    return $s;
}