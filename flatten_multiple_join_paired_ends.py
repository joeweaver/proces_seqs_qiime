#!/usr/bin/python

import argparse
import os
import shutil

__author__ = "Joseph E. Weaver"
__license__ = "MIT"
__version__ = "0.1.1"
__maintainer__ = "Joseph E. Weaver"
__email__ = "jeweave4@ncsu.edu"

################################################################
#Flattens the output of multiple_joined_ends.py from QIIME 1.9.0
################################################################

#create list of files we want to copy and where to copy them
def files_to_copy(read_dir,write_dir):
    sd_pairs=[]
    for root, dirs,files in os.walk(read_dir):
        for file in files:
            if "fastqjoin.join.fastq" == file:
                sd_pairs.append((os.path.abspath(os.path.join(root,file)),
                                os.path.abspath(os.path.join(write_dir,
                                    os.path.basename(os.path.normpath(root))    
                                        +".paired.fastq"))))
    return(sd_pairs)


if __name__ == '__main__':
    #handle command line args
    parser = argparse.ArgumentParser(description='Flatten the directories created by multiple_join_paired_ends.py')
    parser.add_argument("-o","--output_dir",help="directory to output flattened joins",required=True)
    parser.add_argument("-i","--input_dir",help="directory to read",required=True)
    parser.add_argument("-w","--print_only",action='store_true',dest='print_only',help="Print the commands*, but don't execute. *commands are given as if on the commandline, we use python os commands in code",default=False)
    parser.add_argument("-v","--verbose",action='store_true',dest='verbose',help="Be verbose - Right now, just 'X files copied to DIR'",default=False)
    args = parser.parse_args()

    copies_made=0 #track the number of copies made for --verbose
    
    #go through each pair of src and destination files,create the output dir
    #if necessary, then copy the file
    #TODO would be nice if this were atomic
    for src,dst in files_to_copy(args.input_dir,args.output_dir):
        #doing a few gymnastics to avoida  race condition of directory creation        
        if not os.path.exists(args.output_dir):
            try:
                if(args.print_only):
                    print("If it doesn't exist: mkdir "+args.output_dir)
                else:
                    os.mkdir(args.output_dir)
            except OSError as exc:
                if exc.errno != errno.EEXIST:
                    raise
                pass    
        if(args.print_only):
            print("copy "+src+" "+dst)
        else:
            #actually copy the files
            #TODO check for overwrite
            #shouldn't have clobber issues in our toolchain, but wouldn't hurt
            #to be careful. Besides, this gets a lot slower when it tries to
            #copy to an existing file.
            with  open(dst,"w") as f:
                with open(src) as sf:
                    shutil.copyfileobj(sf,f)
                    copies_made=copies_made+1

    #report success if verbose
    if(args.verbose):
        print("Copied " + str(copies_made) + " files to " + args.output_dir)


