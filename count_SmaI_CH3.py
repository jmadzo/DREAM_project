#!/usr/bin/env python
"""
This script asigne bowtie/bwa mapped reads to the SmaI sides and count number of methylated and unmethylated reads
"""
from __future__ import division
__version__ = '0.5'
__author__ = 'hanghang & jmadzo'

import sys
import argparse
import pysam
import os.path
import collections
import subprocess

def extractFileName(input_file_path):
    '''extract file name from shell command when longer path is provided'''
    return os.path.basename(input_file_path)

def splitFileName(file_name):
    '''split file name to name root and extention'''
    return os.path.splitext(file_name)

def isBamSorted(file_name, open_mode):
    '''Opens BAM/Sam file and check if file is sorted '''
    try:
        open_for_sortCheck = pysam.Samfile(file_name, open_mode, check_header=True)
        sortCheck = open_for_sortCheck.header['HD']['SO']
    except Exception, e:
        raise e, "\n Check if file has header"
    finally:
        open_for_sortCheck.close()
        print "Sorting tag:", sortCheck
        return sortCheck

def percent(M,U):
    '''From number of methylated and unmethylated reads calculates precentage of methylC reads handles ZeroDivisionError. No reads output: None'''
    try:
        ratio= M/(M+U)*100
    except ZeroDivisionError:
        ratio=None
    finally:
        return ratio

def main():
    '''main function '''
    ###################### check command line arguments ############################################
    if len(sys.argv)==1: sys.exit('\n no arguments provided, for help: --help')

    ##############################  parsing command line arguments  ################################
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i',action='store', dest='input_file', help='input bam file')
    parser.add_argument('-q', action='store', dest='mapq',type=int, help='filter MAPQ reads quality default=5',default=5)
    parser.add_argument('-v', '--verbose', action='store_true', help='increase output verbosity on the screen')
    
    command_line=parser.parse_args()
    input_file_path=command_line.input_file
    mapq=command_line.mapq
    on_screen=command_line.verbose

    ###### parse file name ######
    file_name = extractFileName(input_file_path)
    file_name_root, file_ext = splitFileName(file_name)

    if file_ext == ".bam": open_mode="rb"
    elif file_ext == ".sam": open_mode="r"
    else: sys.exit("Run stopped, %s is not meaningful file extension" %file_ext)

    ##################################################################################################

    ###### Open file, make sure that bam file is sorted and indexed  ######     
    if isBamSorted(file_name, open_mode)=="unsorted":
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

    ###### when file is already sorted open it, check for index file if not present make one ######
    else :
        if not os.path.isfile(file_name+".bai"):
            print "creating index file\n"
            pysam.index(file_name)
        else: pass
        samfile = pysam.Samfile(file_name, open_mode)

    
    ###### some pysam file stats, this is not important, would be deleted  ######
    #print "samfile.filename:", samfile.filename
    #print "samfile.count():", samfile.count(until_eof=True)
    #print samfile.lengths
    #print samfile.text
    #header = samfile.header
    #print "samfile.header:\n", samfile.header

    print "\nitetare thru reads:"

    ###### open SmaI side file and create ordered dictionary/hash ######
    try:
        f=open("SmaI_sides_keys.txt","r")
        smai_sides=[row.strip().split('\t') for row in f]
        SmaI_DB=collections.OrderedDict([((chromosome,int(position)),({'ID':smai_ID, 'M':0, 'U':0})) for chromosome,position,smai_ID in smai_sides])
    except Exception, e:
        print "Probably no SmaI_sides_keys.txt file in directory "
        raise e
    finally:
        f.close()

    count_reads_low_mapq=0
    count_non_SmalI=0
    count_unmapped = 0
    count_read_pass = 0
    samfile_used_reads = pysam.Samfile("out_used_reads_"+file_name, "wb", template=samfile)
    samfile_not_SmaI_read = pysam.Samfile("out_not_SmaI_read_"+file_name, "wb", template=samfile)
    samfile_reads_low_mapq = pysam.Samfile("out_low_mapq_reads_"+file_name, "wb", template=samfile)
    samfile_unmmaped_reads = pysam.Samfile("out_unmmaped_reads_"+file_name, "wb", template=samfile)
    
    for (counter, read) in enumerate(samfile.fetch(until_eof = True)):
        if counter%1000000==0: print "procesed",counter/1000000,"M reads"
        if read.is_unmapped:
            count_unmapped+=1
            "check_for_spike(read) ... check forspikes"
            samfile_unmmaped_reads.write(read)
            continue

        elif read.is_reverse:
            if read.seq[-3:]=="CCC":
                SmaI_side=read.aend+1
                methylated=False
            elif read.seq[-5:]=="CCCGG":
                SmaI_side=read.aend-1
                methylated=True
            else:
                SmaI_side=None  #save junk -> samfile_not_SmaI_read

        else:
            if read.seq[:3]=="GGG":
                SmaI_side=read.pos+1
                methylated=False
            elif read.seq[:5]=="CCGGG":
                SmaI_side=read.pos+3
                methylated=True
            else:
                SmaI_side=None #save junk -> samfile_not_SmaI_read

        SmaI_side_key=("chr"+str(read.tid),SmaI_side)

        if SmaI_side_key in SmaI_DB:
            if read.mapq > mapq:
                count_read_pass+=1
                if methylated: SmaI_DB[SmaI_side_key]['M']+=1
                if not methylated: SmaI_DB[SmaI_side_key]['U']+=1
                samfile_used_reads.write(read)

            else:
                count_reads_low_mapq+=1
                samfile_reads_low_mapq.write(read)
        else:
            count_non_SmalI+=1
            samfile_not_SmaI_read.write(read) #saving junk file

    samfile.close()
    samfile_used_reads.close()
    samfile_reads_low_mapq.close()
    samfile_not_SmaI_read.close()
    samfile_unmmaped_reads.close()

    print
    count_total = counter+1
    count_mapped = count_read_pass+count_reads_low_mapq+count_non_SmalI # count_total-count_unmapped,
    print count_total, "total number of reads"
    print
    print "{} reads unpapped => {:.3f} %".format(count_unmapped, count_unmapped/count_total*100 )
    print "{} reads mapped => {:.3f} %".format(count_mapped, (count_total-count_unmapped)/count_total*100 )
    print "{} reads above {} MAPQ quality limit => {:.3f} % of mapped or {:.3f} % of total".format(count_read_pass, mapq, count_read_pass/count_mapped*100, count_read_pass/count_total*100 )
    print "%s reads under %s MAPQ quality limit => %.3f %% of mapped  %.3f %% of total" %(count_reads_low_mapq, mapq, count_reads_low_mapq/count_mapped*100, (count_reads_low_mapq/count_total)*100)
    print
    print "%s of mapped reads are NOT in SmaI database => %.4f %% of mapped or %.4f %% of total" %(count_non_SmalI, count_non_SmalI/count_mapped*100 ,(count_non_SmalI/count_total)*100)
    print
    print "count_read_pass + count_reads_low_mapq + count_non_SmalI ->  %s + %s + %s = %s  ==  mapped -> %s" %(count_read_pass, count_reads_low_mapq, count_non_SmalI, (count_read_pass+count_reads_low_mapq+count_non_SmalI) , count_mapped)
    print "mapped + unmapped -> %s + %s = %s == total -> %s" %(count_mapped, count_unmapped , count_mapped+count_unmapped, count_total)
        
    print
    with open('SmaI_sites_output.txt', 'w') as output:
        for key, value in SmaI_DB.iteritems():
            percentage = percent(value['M'], value['U'])
            if on_screen:
                print '{0}\t{1}'.format(*key),"\t{ID}\t{M}\t{U}".format(**value),"\t{}".format(percentage)
            output.write('{1}\t{2}\t{ID}\t{M}\t{U}\t{0}'.format(percentage, *key,**value)+"\n")


#sys.exit("\nOK, good up here\n")
if __name__ == "__main__":
    import timex as tx
    start, startLocal=tx.startTime()
    main()
    tx.stopTime(start,startLocal)