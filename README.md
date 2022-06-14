# RPREP

* RPREP: Ribosomal Profiling Reports, an open-source cloud-enabled framework for reproducible Ribosomal Profiling data processing, analysis, and result reporting.

## INSTALLATION
 
### Option 1 (local Ubuntu server): 

* install LTS Ubuntu desktop version 18.04.02 on your machine (http://releases.ubuntu.com/18.04.02/, install from bootable USB)
* copy/clone the RPREP github source code to your Ubuntu machine (git clone https://github.com/emmesgit/RPREP.git). Ensure read write and execute permissions are set for RPREP (chmod -R u+rwx RPREP).
* execute our installation shell script to install the software on your own Ubuntu machine (sh RPREP/ubuntu/install-software-v2.1.0.sh)

### Option 2 (RPREP AWS AMI):

* initialize RPREP AMI (https://aws.amazon.com, AMI ID: RPREP RSEQREP (Ribosome Profiling and RNA-Seq Reports) v2.1 (ami-00b92f52d763145d3))).  To do this, create a AWS account (https://aws.amazon.com), log into the AWS console, and navigate to the EC2 resources.  Next select AMIs in the navigation pane and select public images.  Finally search for RPREP, find "RPREP RSEQREP (Ribosome Profiling and RNA-Seq Reports) v2.1" and launch the AMI (ensure you are in the US EAST (N. Virginia) region).
* using an ssh command line connection or X2GO GUI supported connection client (https://wiki.x2go.org/doku.php/download:start) (XFCE session type) log into the RPREP AMI using username: repuser and password: repuser2019. The IP address of your instance can be found on the AWS console EC2->Instances->"IPv4 Public IP".  Please note this IP may change every time the machine is started/stopped.
* copy/clone the RPREP github source code to your Ubuntu machine (git clone https://github.com/emmesgit/RPREP.git).  Ensure read write and execute permissions are set for RPREP (chmod -R u+rwx RPREP).

Additional information on AWS and AMI configuration can be found in RPREP/aws/aws_instructions.docx

### Option 3 (RPREP Docker Image):

* On an Ubuntu 18.04.02 machine with Docker software installed (https://docs.docker.com/install), pull RPREP image from the Docker repository (docker pull emmesdock/rseqrep) [https://hub.docker.com/r/emmesdock/rseqrep/]. Run the docker image as a container in interactive mode (docker run --name rseqrep1 -i -t emmesdock/rseqrep /bin/bash).
* copy/clone the RPREP github source code to the container (git clone https://github.com/emmesgit/RPREP.git).  Ensure read write and execute permissions are set for RPREP (chmod -R u+rwx RPREP).

## EXECUTION

* Download GMT formatted gene sets for pathway enrichment (if no gene sets are specified, pathway enrichment is not performed).  We provide a custom script (https://github.com/emmesgit/RPREP/blob/master/source/shell/download-gene-sets.sh) to download Blood Transcription Modules (option btm), Reactome pathways (reactome option), and KEGG pathways (kegg option). Note for the KEGG option you need to either be an academic user or have a commercial KEGG license (http://www.kegg.jp/kegg/rest).  Gene sets can also be downloaded from MSigDB after email registration (http://software.broadinstitute.org/gsea/msigdb).  Specify the locations and labels of the gene sets in the gmt_entrez_files and gmt_entrez_files_labels fields of the worklfow_config tab of the configuration file (RPREP/config/config.xlsx).
* fill out the configuration file (RPREP/config/config.xlsx).  RPREP/case-study/config-henn.xlsx is an example of a complete configuration file.
* execute start-to-end analysis (sh RPREP/run-all.sh --config XXX --threads XXX --log XXX) or for a particular component (sh RPREP/run-pre-processing.sh --config XXX --threads XXX --log XXX; sh RPREP/run-analysis.sh --config XXX --threads XXX --log XXX; sh RPREP/run-report.sh --config XXX --log XXX).
 
## TROUBLESHOOTING

Upon inspecting the initial run of the report, you may find that the configuration option that you initially chose does not fit your data.  For example, you inspect the reverse cumulative distribution function plot comparing log count per million cutoffs with the number of retained genes.  You find that the cutoff you had originally selected resulted in too few genes for the analysis.  To remedy this, you would update the configuration file to reflect a more appropriate log counts per million cutoff.  You then determine the steps of the analysis that will be affected by this configuration change.  In this case, the analysis and report steps are affected.  Re-execute the analysis and report steps (sh RPREP/run-analysis.sh --config XXX --threads XXX --log XXX; sh RPREP/run-report.sh --config XXX --log XXX).  The configuration file is re-parsed each time a RPREP/run-* script is executed.
 
## BUGS/ISSUES

If you identify bugs/issues, please submit them via the GitHub Issue Tracker 
https://github.com/emmesgit/RPREP/issues


## RELEASE NOTES 

### RPREP Ribosomal Profiling Reports - Version 1.1.0

#### New Features:

* Utilizing multicore functionality of cutadapt within the trimadapters rule.  
* Added a --dryrun arg which can passed to the pre-processing execution to print the list of jobs to be completed without actually executing the workflow.

#### Bug Fixes:

* The strandedness and paired end options provided in the configuration file was never properly linked to the feature_counts rule in snakemake.  This functionality is now working properly.

### RPREP Ribosomal Profiling Reports - Version 1.0.1

#### Bug Fixes:

* OpenSSL encryption hash can now be specified as md5 or sha256 in the configuration xlsx file.  previously, md5 was hard-coded.

### RPREP Ribosomal Profiling Reports - Version 1.0.0

* Initial release.

 
