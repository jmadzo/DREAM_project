#!/usr/bin/env python
"""
This script asigne bowtie/bwa mapped reads to the SmaI sides and count number of methylated and unmethylated reads
"""
from __future__ import division
__version__ = '0.1'
__author__ = 'hanghang & jmadzo'

import sys
import argparse
import pysam
import os.path

#extension = os.path.splitext(filename)[1][1:]

def extractFileName(input_file_path):
	'''extract file name from shell command when longer path is provided'''
	return os.path.basename(input_file_path)

def splitFileName(file_name):
	'''split file name to name root and extention'''
	return os.path.splitext(file_name)

def main():
    '''main function '''
    ###################### check command line arguments ############################################
    if len(sys.argv)==1: sys.exit('\n no arguments provided, for help: --help')

    ##############################  parsing command line arguments  ################################
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i',action='store', dest='input_file', help='input bam file')
    parser.add_argument('-q', action="store", dest='mapq',type=int, help='filter MAPQ reads quality',default=10)

    command_line=parser.parse_args()
    input_file_path=command_line.input_file
    mapq=command_line.mapq

    ###### parse file name ######
    file_name = extractFileName(input_file_path)
    file_name_root, file_ext = splitFileName(file_name)

    if file_ext == ".bam": open_mode="rb"
    elif file_ext == ".sam": open_mode="r"
    else: sys.exit("Run stopped, %s is not meaningful file extension" %file_ext)

    ###### first make sure that bam file is sorted and indexed ######
    #check if sorted pass
    open_for_sortCheck = pysam.Samfile(file_name, open_mode, check_header=True)
    sortCheck = open_for_sortCheck.header['HD']['SO']
    open_for_sortCheck.close()
    #print "something went wrong, probably missing header"
    
    print "Sorting tag:", sortCheck, "\n"

    if sortCheck=="unsorted":
    	bam_sorted=file_name_root+"_sort"
    	print " File: %s, \n\t is unsorted, start sorting, this can take a while,\n\n sorted file: %s  \n\t is going to be store in current directory\n" %(file_name, bam_sorted+file_ext)
    	
    	try:
            pysam.sort(file_name, bam_sorted)
            print "creating index file\n"
            pysam.index(bam_sorted+file_ext)

        except Exception as err:
    	    raise err
            pass

        samfile = pysam.Samfile(bam_sorted+file_ext, open_mode)

    ### open the file ###
    else :
    	if not os.path.isfile(file_name+".bai"):
    		print "creating index file\n"
        	pysam.index(file_name)
        else: pass
        samfile = pysam.Samfile(file_name, open_mode)

    
    
    print "samfile.filename:", samfile.filename
    print
    print "samfile.count():", samfile.count('chr21', 10000000, 10200000)
    print
    #print "samfile.header:\n", samfile.header
    
    print
    #sys.exit("\nOK, good up here\n")

    print "\nitetare thru reads\n"
    
    #from itertools import izip

    f=open("SmaI_sides_keys.txt","r")
    smai_sides=[row.strip().split('\t') for row in f]
    f.close()
    Sm_DB={k:v for k,v in smai_sides}

    #Sm_DB=dict(itertools.izip(f.readline().strip().split()))
    #Sm_DB={key: value for key, value in itertools.izip(f.readline().strip().split())}
    
    all_reads={}
    for (counter, read) in enumerate(samfile.fetch('chr21', 10000000, 10200000)):
        all_reads[counter]=read
        if read.is_reverse:
            if read.seq[3]=="GGG":
                SmaI_side=read.aend
                methylated=False
            elif read.seq[5]=="CCGGG":
                SmaI_side=read.aend-2
                methylated=True
            else:
                continue

        else:
            if read.seq[3]=="GGG":
                SmaI_side=read.pos
                methylated=True
            elif read.seq[5]=="CCGGG":
                SmaI_side=read.pos+2
                methylated=True
            else:
                continue
        
        if SmaI_side in SmaI_DB:
            if read.mapq > mapq:


                '''
                print counter
                print "read:",read
                print "chr", read.tid,
                print ":", read.pos,"-",read.aend,
                print "is_reverse:", read.is_reverse
                print "read.seq:\t", read.seq
                print "read.query:\t", read.query
                print "read.alen:",read.alen
        
                print "rnext:", read.rnext,
                print "is_paired:", read.is_paired,
                print "is_proper_pair:", read.is_proper_pair,
                print "isize:", read.isize,
                print "mapq:", read.mapq
                print "\n"
                '''

            else: continue
        else: 
            print read.tid, SmaI_side, "Check check for polymorphism"

    samfile.close()

    sys.exit("\nOK, good up here\n")

    #output_bamfile = pysam.Samfile("out"+imput_bam_sorted, "wb", template=samfile)
    #output_bamfile.close()

if __name__ == "__main__":
	main()