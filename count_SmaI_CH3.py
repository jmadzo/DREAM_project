#!/usr/bin/env python
"""
This script asigne bowtie/bwa mapped reads to the SmaI sides and count number of methylated and unmethylated reads
"""
from __future__ import division
__version__ = '0.5'
__author__ = 'jmadzo'

import sys
import argparse
import pysam
import os.path
import collections

def extractFileName(input_file_path):
    '''extract file name from shell command when longer path is provided'''
    return os.path.basename(input_file_path)

def splitFileName(file_name):
    '''split file name to name root and extention'''
    return os.path.splitext(file_name)

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
    parser.add_argument('-i', '--input_file', action='store', dest='input_file', help='input bam file')
    parser.add_argument('-g','--genome_table',action='store', dest='genome_table', help='table with side for your genome build',default="SmaI_sides_keys_hg19.txt")
    parser.add_argument('-q', '--mapq', action='store', dest='mapq',type=int, help='filter MAPQ reads quality default=5',default=5)
    parser.add_argument('-bed','--ouput_bed', action='store_true', help='ouput table will be in bed like format')
    parser.add_argument('-v', '--verbose', action='store_true', help='increase output verbosity on the screen')
    parser.add_argument('-ob','--ouput_bam', action='store_true', help='ouput BAM files with reads falling into use/low_quqlity/unmap reads' )

    
    command_line=parser.parse_args()
    input_file_path=command_line.input_file
    genome_table=command_line.genome_table
    mapq=command_line.mapq
    bed_like=command_line.ouput_bed
    on_screen=command_line.verbose
    bam_out=command_line.ouput_bam

    ###### parse file name ######
    file_name = extractFileName(input_file_path)
    file_name_root, file_ext = splitFileName(file_name)

    if file_ext == ".bam": open_mode="rb"
    elif file_ext == ".sam": open_mode="r"
    else: sys.exit("Run stopped, %s is not meaningful file extension" %file_ext)

    ##################################################################################################

    ###### Open file ######
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
        f=open(genome_table,"r")
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
    if bam_out:
        samfile_used_reads = pysam.Samfile(file_name+"_out_used_reads", "wb", template=samfile)
        samfile_not_SmaI_read = pysam.Samfile(file_name+"_out_not_SmaI_read", "wb", template=samfile)
        samfile_reads_low_mapq = pysam.Samfile(file_name+"_out_low_mapq_reads", "wb", template=samfile)
        samfile_unmmaped_reads = pysam.Samfile(file_name+"_out_unmmaped_reads", "wb", template=samfile)
    
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
                if bam_out: samfile_used_reads.write(read)

            else:
                count_reads_low_mapq+=1
                if bam_out: samfile_reads_low_mapq.write(read)
        else:
            count_non_SmalI+=1
            if bam_out: samfile_not_SmaI_read.write(read) #saving junk file

    samfile.close()
    try:
        samfile_used_reads.close()
        samfile_reads_low_mapq.close()
        samfile_not_SmaI_read.close()
        samfile_unmmaped_reads.close()
    except: pass

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
    
    ###### Writing stats into the file ###### hack for now but it's ugly, needs to be rewritten
    stats=open(file_name+'stats.txt', 'w')
    stats.write("{} reads unpapped => {:.3f} %".format(count_unmapped, count_unmapped/count_total*100 ))
    stats.write("{} reads mapped => {:.3f} %".format(count_mapped, (count_total-count_unmapped)/count_total*100 )
    stats.write("{} reads above {} MAPQ quality limit => {:.3f} % of mapped or {:.3f} % of total".format(count_read_pass, mapq, count_read_pass/count_mapped*100, count_read_pass/count_total*100 )
    stats.write("%s reads under %s MAPQ quality limit => %.3f %% of mapped  %.3f %% of total" %(count_reads_low_mapq, mapq, count_reads_low_mapq/count_mapped*100, (count_reads_low_mapq/count_total)*100)
    stats.write("")
    stats.write("%s of mapped reads are NOT in SmaI database => %.4f %% of mapped or %.4f %% of total" %(count_non_SmalI, count_non_SmalI/count_mapped*100 ,(count_non_SmalI/count_total)*100)
    stats.write("")   
    stats.write("count_read_pass + count_reads_low_mapq + count_non_SmalI ->  %s + %s + %s = %s  ==  mapped -> %s" %(count_read_pass, count_reads_low_mapq, count_non_SmalI, (count_read_pass+count_reads_low_mapq+count_non_SmalI) , count_mapped)
    stats.write("mapped + unmapped -> %s + %s = %s == total -> %s" %(count_mapped, count_unmapped , count_mapped+count_unmapped, count_total)
    stats.close()


    print
    with open(file_name+'_SmaI_sites.txt', 'w') as output:
        for key, value in SmaI_DB.iteritems():
            percentage = percent(value['M'], value['U'])
            if on_screen:
                print '{0}\t{1}'.format(*key),"\t{ID}\t{M}\t{U}".format(**value),"\t{}".format(percentage)
            if bed_like:
                output.write('{1}\t{2}\t{ID}\t{M}\t{U}\t{0}'.format(percentage, *key,**value)+"\n")
            else:
                output.write('{ID}\t{1}\t{2}\t{M}\t{U}\t{0}'.format(percentage, *key,**value)+"\n")


#sys.exit("\nOK, good up here\n")
if __name__ == "__main__":
    import timex as tx
    start, startLocal=tx.startTime()
    main()
    tx.stopTime(start,startLocal)