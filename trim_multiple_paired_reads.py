#!/usr/bin/python

import argparse
import os
import shutil
import subprocess

__author__ = "Joseph E. Weaver"
__license__ = "MIT"
__version__ = "0.1.1"
__maintainer__ = "Joseph E. Weaver"
__email__ = "jeweave4@ncsu.edu"

################################################################
#Calls Trimmomatic on multiple files
################################################################

#borrows heavily from flatten_multiple_paired_ends
#TODO refactor out similarities?

#create list of files we want to trim and where to write them
#in our current pipeline, this assumes we're looking at the mulitdirectory
#output from multiple_extract_barcodes.py
def files_to_trim(read_dir,write_dir):
    sd_pairs=[]
    for root, dirs,files in os.walk(read_dir):
        for file in files:
            if "reads.fastq" == file:
                sd_pairs.append((os.path.abspath(os.path.join(root,file)),
                                os.path.abspath(os.path.join(write_dir,
                                    os.path.basename(os.path.normpath(root))    
                                        +".trm_dbc_paired.fastq"))))
    return(sd_pairs)


if __name__ == '__main__':
    #handle command line args
    parser = argparse.ArgumentParser(description='Flatten the directories created by multiple_join_paired_ends.py')
    parser.add_argument("-o","--output_dir",help="directory to output trimmed reads",required=True)
    parser.add_argument("-i","--input_dir",help="directory to read",required=True)
    parser.add_argument("-w","--print_only",action='store_true',dest='print_only',help="Print the commands, but don't execute.",default=False)
    parser.add_argument("-v","--verbose",action='store_true',dest='verbose',help="Be verbose - Right now, just 'X files copied to DIR'",default=False)
    args = parser.parse_args()

    files_trimmed=0 #track the number of files_trimmed for --verbose
    
    #create a list of commands invoking trimmomatic on each file
    commands=[]
    for src,dst in files_to_trim(args.input_dir,args.output_dir):
        #unlike flatten_multiple_paired_ends.py, we assume that the directory has been
        #created by the pipleine just prior to calling this.       


        #create the list of calls to trimmotatic
        #we assume it exsits and is callable
        #~/bin should be PATH at this point, be we won't assume it
        #TODO handle errors from failed calls/runs and file clobbers
        command = ("java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE",
                   src,
                   dst,
                   "LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36",
					" >> /dev/null")
        commands.append(command)
        files_trimmed=files_trimmed+1
    for command in commands:
        if(args.print_only):
            print(' '.join(command))                
        else:
            #sending this as a single string is somewhat insecure
            #however, we're not taking input from untrusted users
			#sending output to dev/null
            FNULL = open(os.devnull,'w')
            subprocess.call(' '.join(command),shell=True,stdout=FNULL,stderr=subprocess.STDOUT)
    #report success if verbose
    if(args.verbose):
        print("Trimmed " + str(files_trimmed) + " files to " + args.output_dir)


