#!/usr/bin/env python
"""
This script assigns bowtie/bwa mapped reads from BAM/SAM file to the all SmaI sides and count number of methylated and unmethylated reads.

standard ouput format:
SmaI_ID     chrN    pos     M   U   ratio
SmaI_hg19_1 chr1    10498   1   0   100.0
SmaI_hg19_2 chr1    20207   0   0   23.789
SmaI_hg19_3 chr1    24074   0   0   None

bed format ouput:
chrN    star    end     ratio   SmaI_ID     M   U  
chr1    10496   10498   100.0   SmaI_hg19_1 1   0
chr1    20205   20207   None    SmaI_hg19_2 0   0
"""
from __future__ import division
__version__ = '1.0'
__author__ = "JP Issa's lab"

import sys
import argparse
import pysam
import os.path
import collections
import re
from subprocess import Popen

def extractFileName(input_file_path):
    '''extract file name from shell command when longer path is provided'''
    return os.path.basename(input_file_path)

def splitFileName(file_name):
    '''split file name to name root and extention'''
    return os.path.splitext(file_name)

def percent(M,U):
    '''From number of methylated and unmethylated reads calculates precentage of methylC reads handles ZeroDivisionError. No reads output: None'''
    try:
        ratio= round(M/(M+U)*100,2)
    except ZeroDivisionError:
        ratio=None
    finally:
        return ratio

def spike_checker(read, samfile_spike_reads, bam_out=False):
    '''Check if unmapped read sequence is in the spike database, if yes increase the methylation counter of match spike'''
    read_sequence=read.seq
    if read_sequence[:5]=="CCGGG":
        status="M"
        sequence_query=read_sequence[2:]
    elif read_sequence[:3]=="GGG":
        status="U"
        sequence_query=read_sequence
    else:
        status=None
        return None
    
    seq=re.compile(sequence_query)
    for spike in db_spikes:
        if seq.search(db_spikes[spike].get('seq')):
            db_spikes[spike][status]+=1
            if bam_out: samfile_spike_reads.write(read)

    return None

def main():
    '''main function '''
    ###################### check command line arguments ############################################
    if len(sys.argv)==1: sys.exit(__doc__+'\n no arguments provided, for help: --help')

    ##############################  parsing command line arguments  ################################
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i', '--input_file', action='store', dest='input_file', help='input bam file')
    parser.add_argument('-g','--genome',action='store', dest='genome_table', help='table with RE side position for accessed genome build',default="SmaI_sites_keys_hg19.txt")
    parser.add_argument('-q', '--mapq', action='store', dest='mapq',type=int, help='filter MAPQ reads quality default=5',default=5)
    parser.add_argument('-bed','--ouput_bed', action='store_true', help='ouput table will be in bed like format')
    parser.add_argument('-v', '--verbose', action='store_true', help='increase output verbosity on the screen')
    parser.add_argument('-bam','--ouput_bam', action='store_true', help='ouput BAM files with reads falling into use/low_quality/unmap reads' )
    parser.add_argument('-s','--silent', action='store_false', help='run in silent mode (no screen output)')
    parser.add_argument('-n', '--no_header', action='store_false', help='prints file output without header')
    parser.add_argument('-spike', action='store', dest='spike_file', help='imput file with spike name and sequencies', default=None)
    parser.add_argument('-bowtie', action='store', dest='bowtie_par', help='add bowtie parameters input files and output file include=')

    command_line=parser.parse_args()
    input_file_path=command_line.input_file
    genome_table=command_line.genome_table
    mapq=command_line.mapq
    bed_like=command_line.ouput_bed
    on_screen=command_line.verbose
    bam_out=command_line.ouput_bam
    silent=command_line.silent
    no_header=command_line.no_header
    spike_file = command_line.spike_file
    bowtie_parameters=command_line.bowtie_par

    header=(True and no_header)
    print bool(spike_file)

    ###### Bowtie call ########
    if bowtie_parameters:
        command_bowtie= "bowtie2 " + bowtie_parameters #+ " | samtools view -bS - | samtools sort - fighetta.bam"
        print command_bowtie
        bowie_output=re.compile('\w+.[bs]am')
        input_file_path = re.search(bowie_output,bowtie_parameters).group()
        print "\nexpected output file: %s\n" % input_file_path

        try:
            proc = Popen(command_bowtie, shell=True)
            proc.wait()
            print "done"
            print
            
        except Exception, e:
            print "bowtie did't run"
            raise e
            sys.exit("\nbowtie did't run\n")

    ###### parse file name ######
    file_name = extractFileName(input_file_path)
    file_name_root, file_ext = splitFileName(file_name)

    if file_ext == ".bam": open_mode="rb"
    elif file_ext == ".sam": open_mode="r"
    else: sys.exit("Run stopped, %s is not meaningful file extension" %file_ext)

    ##################################################################################################

    ###### Open file ######
    try:
        samfile = pysam.Samfile(file_name, open_mode)
    except IOError as ioe:
        sys.exit("BAM/SAM file error: {}\nCheck spelling of file name \nIf you ran -bowtie option, there is good chance it didn't run".format(ioe))
    
    ###### some pysam file stats, this is not important, would be deleted  ######
    #print "samfile.filename:", samfile.filename
    #print "samfile.count():", samfile.count(until_eof=True)
    #print samfile.lengths
    #print samfile.text
    #header = samfile.header
    #print "samfile.header:\n", samfile.header

    if silent: print "\nitetare thru reads:"

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

    ###### open spikes file and create spike DB ######
    global db_spikes
    try:
        spk=open(spike_file,"r")
        spike_list=[line.strip().split() for line in spk]
        db_spikes={name:{'seq':sequence,'M':0,'U':0} for name,sequence in spike_list}
    except:
        db_spikes = {'None': {'M': 0,'U': 0,'seq': 'NNNNNNNNNNNNNNN..........NNNNNNNNNNNNNNN'}}

    if bam_out:
        day_stapm="_"+tx.today()
        path_DIR=file_name_root+"_BAM_"+day_stapm+"/"
        if not os.path.exists(path_DIR): os.mkdir(path_DIR)
        samfile_used_reads = pysam.Samfile(path_DIR+file_name_root+"_out_used_reads.bam", "wb", template=samfile)
        samfile_not_SmaI_read = pysam.Samfile(path_DIR+file_name_root+"_out_not_SmaI_read.bam", "wb", template=samfile)
        samfile_reads_low_mapq = pysam.Samfile(path_DIR+file_name_root+"_out_low_mapq_reads.bam", "wb", template=samfile)
        samfile_unmmaped_reads = pysam.Samfile(path_DIR+file_name_root+"_out_unmmaped_reads.bam", "wb", template=samfile)
    if bool(spike_file) and bam_out:
        samfile_spike_reads = pysam.Samfile(path_DIR+file_name_root+"_out_spike_reads.bam", "wb", template=samfile)
    else:
        samfile_spike_reads = None
    
    for (counter, read) in enumerate(samfile.fetch(until_eof = True)):
        if silent and counter%1000000==0: print "processed",counter/1000000,"M reads"
        if read.is_unmapped:
            count_unmapped+=1
            ''' spike_checker function check if read is in spike DB, if yes increase spike's M or U counter'''
            if bool(spike_file): spike_checker(read, samfile_spike_reads, bam_out)
            if bam_out: samfile_unmmaped_reads.write(read)
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

        SmaI_side_key=(samfile.getrname(read.tid),SmaI_side)

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
    
    #print "last read key: {}\t{}".format(*SmaI_side_key)
    #print "read:", read
    #print "read.tid:", read.tid
    #print "read.rlen:", read.rlen
    #print "read.rname:", read.rname
    #print "read.rnext:", read.rnext
    #print "read.mrnm:", read.mrnm
    #print "samfile.getrname(read.tid):", samfile.getrname(read.tid)
    
    samfile.close()
    try:
        samfile_used_reads.close()
        samfile_reads_low_mapq.close()
        samfile_not_SmaI_read.close()
        samfile_unmmaped_reads.close()
        samfile_spike_reads.close()
    except: pass

    if silent: print
    count_total = counter+1
    count_mapped = count_read_pass+count_reads_low_mapq+count_non_SmalI # count_total-count_unmapped,

    if silent: 
        print "{:,} total number of reads".format(count_total)
        print
        print "{:,} reads unpapped => {:.3f} %".format(count_unmapped, count_unmapped/count_total*100 )
        print "{:,} reads mapped => {:.3f} %".format(count_mapped, (count_total-count_unmapped)/count_total*100 )
        print "{:,} reads above {} MAPQ quality limit => {:.3f} % of mapped or {:.3f} % of total".format(count_read_pass, mapq, count_read_pass/count_mapped*100, count_read_pass/count_total*100 )
        print "{:,} reads under {} MAPQ quality limit => {:.3f} % of mapped  {:.3f} % of total".format(count_reads_low_mapq, mapq, count_reads_low_mapq/count_mapped*100, (count_reads_low_mapq/count_total)*100)
        print "{:,} of mapped reads are NOT in SmaI database => {:.4f} % of mapped or {:.4f} % of total".format(count_non_SmalI, count_non_SmalI/count_mapped*100 ,(count_non_SmalI/count_total)*100)
        print
        print "reads pass mapq:\t{} \nreads low mapq:\t{} \nnon SmaI:\t{} \n".format(count_read_pass, count_reads_low_mapq, count_non_SmalI)
        print "mapped\t{} \nunmapped:\t{} \ntotal:\t{} \n".format(count_mapped, count_unmapped, count_total)
        print
    
    ###### Writing stats into the file ###### hack for now but it's ugly, needs to be rewritten
    stats=open(file_name_root+'_stats.txt', 'w')
    stats.write("{} reads unpapped => {:.3f} %\n".format(count_unmapped, count_unmapped/count_total*100 ))
    stats.write("{} reads mapped => {:.3f} %\n".format(count_mapped, (count_total-count_unmapped)/count_total*100 ))
    stats.write("{} reads above {} MAPQ quality limit => {:.3f} % of mapped or {:.3f} % of total\n".format(count_read_pass, mapq, count_read_pass/count_mapped*100, count_read_pass/count_total*100 ))
    stats.write("{} reads under {} MAPQ quality limit => {:.3f} % of mapped  {:.3f} % of total\n".format(count_reads_low_mapq, mapq, count_reads_low_mapq/count_mapped*100, (count_reads_low_mapq/count_total)*100))
    stats.write("{} of mapped reads are NOT in SmaI database => {:.4f} % of mapped or {:.4f} % of total\n".format(count_non_SmalI, count_non_SmalI/count_mapped*100 ,(count_non_SmalI/count_total)*100))
    stats.write("\n")   
    stats.write("reads pass mapq:\t{:} \nreads low mapq:\t{:} \nnon SmaI:\t{:} \n".format(count_read_pass, count_reads_low_mapq, count_non_SmalI))
    stats.write("mapped:\t{:} \nunmapped:\t{:} \ntotal:\t{:} \n".format(count_mapped, count_unmapped, count_total))
    stats.close()

    with open(file_name_root+'_SmaI_sites.txt', 'w') as output:
        if header and bed_like:
            output.write('chr\tstart\tendt\percent\tSmaI_ID\t5mC\tC\t{4}\ttotal\n')
        elif header:
            output.write('SmaI_ID\tchr\tpos\t5mC\tC\tpercent\ttotal\n')
        else: pass

        for key, value in SmaI_DB.iteritems():
            percentage = percent(value['M'], value['U'])
            coverage = value['M']+value['U']
            if bed_like:
                chrom,pos = key
                output.write('{0}\t{1}\t{2}\t{3}\t{ID}\t{M}\t{U}\t{4}\n'.format(chrom, pos-2, pos, percentage, coverage,**value)) # now prints CpG to access CCCGGG: pos-4, pos+2
            else:
                output.write('{ID}\t{2}\t{3}\t{M}\t{U}\t{0}\t{1}\n'.format(percentage, coverage, *key,**value))
            if on_screen:
                print '{0}\t{1}'.format(*key),"\t{ID}\t{M}\t{U}".format(**value),"\t{}\t{}".format(percentage, coverage)
    if bool(spike_file):
        with open(file_name_root+'_SmaI_spikes.txt', 'w') as spike_output:
            if header:
                spike_output.write('{0:10}\t{1:6}\t{2:6}\t{3:6}\n'.format('Spike_Name ','5mC','C','%'))
            if silent:
                print '{0:10}\t{1:6}\t{2:6}\t{3:6}'.format('Spike_Name','5mC','C','%')

            for key,value in db_spikes.items():
                percentage = percent(value['M'], value['U'])
                spike_output.write('{0:10}\t{M:6d}\t{U:6d}\t{1}\n'.format(key, percentage, **value))
                if silent:
                    print '{0:10}\t{M:6d}\t{U:6d}\t{1}'.format(key, percentage, **value)

#sys.exit("\nOK, good up here\n")
if __name__ == "__main__":
    import timex as tx
    #start, startLocal=tx.startTime()
    main()
    # tx.stopTime(start,startLocal)