#!/usr/bin/python3
"""
This script will use the following input files:
1. input SNP bed file
2. SNP database url (can be found at: http://hgdownload.soe.ucsc.edu/gbdb/hg38/snp/)
3. file that keeps track of which sequences are used "fasta_input.txt"
4. file that keeps track of which sequences are used on both ends "fatatrack_input.txt"
This input will be used to look up the common SNPs in the region of the fusionRNA.
The output will be a list of positions that should be avoided in the template sequence.

Output file:
1. output SNP bed file with positions to avoid
"""
####################################################################################################
####################################  1. Parse arguments  ##########################################
####################################################################################################

# # import all libraries
import os
import argparse

parser = argparse.ArgumentParser(description='give arguments to main snp script')
parser.add_argument('-i', nargs=1, required=True, help='input SNP bed file, example:input_SNP_fusion0.bed, 0-based')
parser.add_argument('-u', nargs=1, required=True, help='SNP url', metavar='SNP url') # when using a differente species than human, the correct SNP database url should be provided; alternatively, this paramater can be set to 'off' if no SNP database is available
parser.add_argument('-f', nargs=1, required=True, help='file that keeps track of which sequences are used "fasta_input.txt"')
parser.add_argument('-t', nargs=1, required=True, help='file that keeps track of which sequences are used on both ends "fatatrack_input.txt"')

args = parser.parse_args()
input_bed = args.i[0]
SNP_url = args.u[0]

####################################################################################################
###################################  2. process the file  #########################################
####################################################################################################

input_file = open(input_bed)
fusionRNA = input_file.read()

# retreive chrom start end info and add nr of bases bases
fusion_chrom1 = fusionRNA.split()[0]
fusion_start1 = int(fusionRNA.split()[1]) # 0-based
fusion_end1 = int(fusionRNA.split()[2]) # 0-based
fusion_chrom2 = fusionRNA.split()[3]
fusion_start2 = int(fusionRNA.split()[4]) # 0-based
fusion_end2 = int(fusionRNA.split()[5]) # 0-based
fusion_ID = fusionRNA.split()[6]
length = int(fusionRNA.split()[7])

# create an output file
output_file = open("output_SNP_" + fusion_ID + ".bed", "w")

# add the arguments for the SNP search in each side of the fusionRNA
arguments = [(fusion_chrom1,fusion_start1,fusion_end1,"left"), (fusion_chrom2,fusion_start2,fusion_end2,"right")]

####################################################################################################
###################################     3. SNP_url "off"   #########################################
####################################################################################################

# Ignore common SNPs in the region of the fusionRNA
if SNP_url == 'off':
	avoid_range = '[]'

####################################################################################################
###################################     3. SNP_url "url"   #########################################
####################################################################################################

# with the url to the database SNPs are searched in the left and right region (looped twice)

else:
    SNP_all = []
    avoid_range = []
    # loop once for the left side and once for the right side
    for arg in arguments:
        chromosome = arg[0]
        start = arg[1]
        end = arg[2]
        position = arg[3]
        # retreive SNPs in that region
        os.system("bigBedToBed " + SNP_url + " -chrom="+ chromosome + " -start=" + str(start) + " -end=" + str(end) + " " + fusion_ID + "_"+ position+ "_tmp.bed")

        # get info on which sequences were used
        fasta = open(args.f[0])
        fasta_track = open(args.t[0])

        # possible SNPs in left and right side
        poss_end = []
        poss_start = []

        # get the positions of the fusionRNA in the template sequence
        for seq, seq_type in zip(fasta, fasta_track):
            
            chrom, start, end = seq.rstrip().replace(":", "-").split("-")
            start = int(start) - 1  	# 1-based to 0-based
            end = int(end) -1           # 1-based to 0-based

            if seq_type[0:4] == '>end':
                poss_end = list(range(start,end)) + poss_end

            elif seq_type[0:6] == '>start':
                poss_start = poss_start + list(range(start,end))

        # this variable is a list of all fusionRNA positions in the right order for template seq
        poss = poss_end + poss_start	
        

        # make list of SNPs in that region
        SNPs = open(fusion_ID + "_"+position+"_tmp.bed")
        
        # get all SNPs in the region
        for SNP in SNPs:
            SNP_pos = int(SNP.split()[1])
            SNP_length = int(SNP.split()[2]) - int(SNP.split()[1])
        
            for pos in range(SNP_length):
                SNP_all.append(SNP_pos)
                SNP_pos += 1




    # retrieve those SNPs that are actually in the template sequence
    
    for snp in SNP_all:
        if snp in poss:
            avoid_range.append(poss.index(snp))  # index can be used as positions in poss are in correct order (they are 0-based and index as well so no correction needed)

# write the positions to avoid in the output file
output_file.write(str(avoid_range))

# close all files
output_file.close()
input_file.close()
SNPs.close()
fasta.close()
fasta_track.close()
