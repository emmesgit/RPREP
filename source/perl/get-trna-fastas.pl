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
# Program:  get-trna-fastas.pl
# Version:  RPREP 1.0.1
# Author:   Travis L. Jensen, Johannes B. Goll, Sami R. Cherikh, and William F. Hooper
# Purpose:  get Ensembl tRNAs
# Input:    N/A
# Output:   ARGV[0]
#           dirname(ARGV[0])/ensembl-trnas-version.txt
#############################################################################################################

use strict;
use warnings;
use File::Basename;
use Time::Piece;

use lib "$ARGV[2]";
use lib "$ARGV[1]/ensembl/modules";
use lib "$ARGV[1]/ensembl-compara/modules";
use lib "$ARGV[1]/ensembl-variation/modules";
use lib "$ARGV[1]/ensembl-funcgen/modules";
use lib "$ARGV[1]/ensembl-io/modules";


use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::ApiVersion;
use Bio::EnsEMBL::DBSQL::SliceAdaptor;
use Bio::EnsEMBL::Slice;
my $ens_version = software_version();
print "[trna download] Using Ensembl Version $ens_version\n";

## argv[0] output fasta location
my $outfile = $ARGV[0];
my $logfile = join "", dirname($outfile), "/ensembl-trnas-version.txt";


## Establish connection to the database, grab adaptor
my $registry = 'Bio::EnsEMBL::Registry';


$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org', # options: 'useastdb.ensembl.org', 'ensembldb.ensembl.org'
    -user => 'anonymous',
);

my $slice_adaptor = $registry->get_adaptor( 'Human', 'Core', 'Slice');


## Open logfile for writing, add today's date and the ensembl version used
my $date = localtime->strftime('%m/%d/%Y');
open(my $lg, '>', $logfile) or die "Could not open $logfile: $!"; 
print $lg "Ensembl version $ens_version  \nDate downloaded: $date \n";

## Loop through chromosomes and export tRNA sequences to fasta file
open(my $f, '>', $outfile) or die "Could not open $outfile: $!"; 
for(my $chromosome=1; $chromosome<23; $chromosome++){
    my $slice = $slice_adaptor->fetch_by_region( 'chromosome', "$chromosome" );
    foreach my $simple_feature ( @{ $slice->get_all_SimpleFeatures('tRNAscan') } ) {

        my $id = $simple_feature->display_id();
        my $seq = $simple_feature->seq();
        print $f ">$id\n";
        print $f $seq."\n"
    }
}

close $f;
close $lg;
