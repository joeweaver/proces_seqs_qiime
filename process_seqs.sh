#!/bin/sh

#####################
# Amplicon metagenomics pipeline using QIIME 1.9.0
# Written by Joe Weaver (jeweave4@ncsu.edu) for the 
# de los Reyes lab #TODO dlr website
# Version 0.1.1
# Latest version at github ()
# Requires helper scripts from 'qiime_helpers' ()
# Based on the manual workflow developed by Ling Wang (lwang18@ncsu.edu)
#
# License Apache 2.0
# Other labs, feel free to use at your own risk
#
# Starts with a tar.gz of MiSeq sequences from our GSL (should be just the unpaired fasta files)
# and ends up with an OTU table, phylo tree, and the original mapping file in /results 
#
#
# Usage:
# Create a new directory containing:
# - mapping.txt <- QIIME metadata mapping file
# -- This should contain the metadata for the samples you will be processing.
# - extract_barcodes.params.txt <- Barcode info for trimming
# -- IMPORTANT: be sure to update the lengths if you don't use our usual primers
# - split_libraries_fastq.params.txt <- params to split_libraries
# -- you do not need to edit this
#
# You can use 'new_analysis.sh <analysis_directory> to create a new directory containing templates for those files
#
# Examples of all of the above files are locally in ~/example_qiime_files/
# They are also in /example_files/ in the git repo
#
# Your sequences (the tar.gz file you got from the GSL) should live somewhere
# else.  For example, my test file 'FSS_small.tar.gz' lives in ~/joe/FSS_small.tar.gz'
#
# Run the script:  process_seqs.sh <seqs.tar.gz>
# 
# For example:
# '~/bin/process_seqs.sh ~/joe/FSS_small.tar.gz'
#
# The mapping file, phylo tree, and biom file will be copied to ./results
# 
# This will likely take a while, you can use 'htop' in a new command window
# to see if things are running.
#####################

#####################
#First some staging
#####################

#Check for our prerequisites
# 1. mapping file
# 2. params files
# 3. path to  sequences

/bin/echo -e "\e[1mChecking prerequisite files and directories\e[0m"

#check for the existence of the mapping file
#TODO arg to specify mapping file?
if [ ! -e mapping.txt ]; then
	/bin/echo -e '\e[31mI could not find \e[33m"mapping.txt"\e[31m in \e[33m'  $PWD '\e[0m'
	exit 10
fi

#validate the mapping file and exit if there are problems
map_check=`validate_mapping_file.py -m mapping.txt -o debug_mapping | awk '{print $1}'`
# First word from validate_mapping_file.py for bad mapping is "Errors" so we're basing it on that
if [ $map_check = "Errors" ]; then
	/bin/echo -e '\e[31mMapping file has errors. Check \e[33m./debug_mapping\e[31m for hints\e[0m'
	exit 11
fi

#check for the existence of params files
#TODO put these in a list and loop thru?
if [ ! -e extract_barcodes.params.txt ]; then
	/bin/echo -e '\e[31mI could not find \e[33m"extract_barcodes.params.txt"\e[31m in \e[33m' $PWD '\e[0m'
	exit 12
fi

#TODO put these in a list and loop thru?
if [ ! -e split_libraries_fastq.params.txt ]; then
	/bin/echo -e '\e[31mI could not find \e[33m"split_libraries_fastq.params.txt"\e[31m in \e[33m' $PWD '\e[0m'
	exit 12
fi

#TODO nicer error message if seqs file not found. Help page?
if [ ! -e $1 ]; then
	/bin/echo -e '\e[31mMissing sequence file \e[33m' $PWD '\e[0m'
	exit 13
fi

#Pre-reqs met, let's make directories and unpack the sequences.

#Create directories for output, if not created by the scripts we'll call
/bin/echo -e '   Creating directories'
mkdir -p '0_raw_fastq'
mkdir -p '3_trimmed_debarcoded_paired'

/bin/echo -e '   Extracting sequences from \e[33m' $1 '\e[0m \n   This may take a moment.'
#arch_size = `du -sk tar $1 | cut -f 1`

tar -xzf $1 -C ./0_raw_fastq

#########################################################################################################
#We're all set to go.  First, take our sequences and prep them for actual analysis
#########################################################################################################
/bin/echo -e "\e[1mPrepping sequences\e[0m"

#TODO allow option for single end workflows? maybe run all?

#joining paired ends.
/bin/echo -e '   Joining paired ends'
multiple_join_paired_ends.py -i 0_raw_fastq -o 1_paired --read1_indicator '_R1' --read2_indicator '_R2'
if [ ! $? -eq 0 ]; then
	exit 13
fi

#basic QIIME file creates a set of directories each with filenames.  This help script flattens it.
fmjpe=`~/bin/flatten_multiple_join_paired_ends.py -i 1_paired -o 1_flattened_paired -v`
/bin/echo -e '  ' $fmjpe

#remove the barcodes and trim
/bin/echo -e '   Removing barcodes and trimming'
multiple_extract_barcodes.py -i 1_flattened_paired -o 2_debarcoded_paired -p extract_barcodes.params.txt
#TODO get rid of some of the stuff trimmomatic outputs
tmpr=`~/bin/trim_multiple_paired_reads.py -i 2_debarcoded_paired -o 3_trimmed_debarcoded_paired -v`
/bin/echo -e '  ' $tmpr

/bin/echo -e '   Merging processed sequences into single file'
multiple_split_libraries_fastq.py -i 3_trimmed_debarcoded_paired -o 4_single_fasta -p split_libraries_fastq.params.txt

#########################################################################################################
#Now that we have our sequences set up, lets start picking OTUS, chimera checking, etc.
#########################################################################################################
/bin/echo -e "\e[1mRunning analysis\e[0m"

#picking open reference OTUS
#TODO the escape codes aren't working
/bin/echo -e '   Picking OTUS using the \e[33mSILVA\e[0m reference DB \e[31m\e[1mThis could take a while...\e[0m' 
#TODO could allow closed reference or complete de novo, but I'd like to see justificiation for it

#rule of thumb is each job (-0#) uses about 4 gigs. -O 9 has been chosen
#based on what's avail after the VM 'tax' and on historical performance

#suppresing taxonomy since we're doing that later
#TODO allow switching reference db?
pick_open_reference_otus.py -i 4_single_fasta/seqs.fna -o 5_uclust_open_ref_otus -s 0.1 -a -O 9 --suppress_taxonomy_assignment -r ~/refdbs/SILVA_128_QIIME_release/rep_set/rep_set_16S_only/97/97_otus_16S.fasta


#assigning taxonomy
/bin/echo -e '   Assigning taxonomy using \e[33mMiDAS\e[0m reference DB. \e[31m\e[1mThis could take a while...\e[0m'
#TODO choice of db as param/parallel DB use

parallel_assign_taxonomy_rdp.py -c 0.8 -i 5_uclust_open_ref_otus/rep_set.fna -o 5_uclust_open_ref_otus/midas_assigned_taxonomy -r ~/refdbs/midas_s123_213/MiDAS_S123_2.1.3.fasta -t ~/refdbs/midas_s123_213/MiDAS_S123_2.1.3.tax -O 4 --rdp_max_memory 10000 --rdp_classifier_fp /qiime_software/rdpclassifier-2.2-release/rdp_classifier-2.2.jar

###################
#Begin chimera removal
###################

/bin/echo -e '   Identifying putative chimeras. \e[31m\e[1mThis could take a while...\e[0m'
#using silva for the aligned reference seqs. MiDAS does not have a Pynast aligned core, but it is a Silva deriviative
#TODO we could use parallel_align_seqs on the midas refdb and rerun chimera slayer to see what's up
parallel_identify_chimeric_seqs.py -i 5_uclust_open_ref_otus/pynast_aligned_seqs/rep_set_aligned.fasta -o 5_uclust_open_ref_otus/pynast_aligned_seqs/chimeric_seqs.txt -m ChimeraSlayer -a ~/refdbs/SILVA_128_QIIME_release/core_alignment/core_alignment_SILVA128.fna -O 8

/bin/echo -e '   Filtering out putative chimera sequences'
#TODO look into this being the correct dir for the output to live in
filter_fasta.py -f 5_uclust_open_ref_otus/pynast_aligned_seqs/rep_set_aligned.fasta -o 5_uclust_open_ref_otus/pynast_aligned_seqs/rep_set_aligned_no_chmc.fasta -s 5_uclust_open_ref_otus/pynast_aligned_seqs/chimeric_seqs.txt -n

filter_alignment.py -i 5_uclust_open_ref_otus/pynast_aligned_seqs/rep_set_aligned_no_chmc.fasta -o 5_uclust_open_ref_otus/pynast_aligned_seqs/

###################
#end chimera removal
###################

##########################
#ready to make out desired output - a tree and biom files
##########################
echo '   Creating phylo tree'
make_phylogeny.py -i 5_uclust_open_ref_otus/pynast_aligned_seqs/rep_set_aligned_no_chmc_pfiltered.fasta -o 5_uclust_open_ref_otus/rep_set_no_chmc.tre

echo '   Creating BIOM file'
#TODO normalize stuff so output is the last arg when possible
make_otu_table.py -i 5_uclust_open_ref_otus/final_otu_map_mc2.txt -o 5_uclust_open_ref_otus/otu_table_w_rdp_no_chmc.biom -t 5_uclust_open_ref_otus/midas_assigned_taxonomy/rep_set_tax_assignments.txt -e 5_uclust_open_ref_otus/pynast_aligned_seqs/chimeric_seqs.txt

mkdir -p 'results'
cp ./5_uclust_open_ref_otus/rep_set_no_chmc.tre ./results/
cp ./5_uclust_open_ref_otus/otu_table_w_rdp_no_chmc.biom ./results/
cp mapping.txt ./results/

/bin/echo -e '\e[32m\e[1mFinished! Mapping file, tree, and biom file are in results\e[0m'

if [ $2 = "e_egg" ]; then 

e_egg="CuKZqyBJIGp1c3QgbG92ZSBzY2FubmluZyBmb3IgbGlmZWZvcm1zLiBMaWZlZm9ybXMuLi4gCuKZqyB5b3UgdGlueSBsaXR0bGUgbGlmZWZvcm1zLi4uIArimasgeW91IHByZWNpb3VzIGxpdHRsZSBsaWZlZm9ybXMuLi4gCuKZqyB3aGVyZSBhcmUgeW91Pwo="

/bin/echo -e '\e[93m\e[5m' 
/bin/echo $e_egg | base64 -d
/bin/echo -e '\e[0m'
#/bin/echo -e '\e[5m\e93m' $msg '\e[0m'

fi

#TODO write snippet for citing software & tracking version/db used
#TODO write snippet for phyloseq getting started
exit 0

#TODO - smarter error checking at each stage?
#TODO - notify when done? - write token to windows shared folder, have script on windows-side tweet/email/message?
#TODO don't redo work if it already has been done (a la make, luigi, or some other workflow script)

