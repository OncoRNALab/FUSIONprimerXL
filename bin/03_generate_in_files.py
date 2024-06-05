#!/usr/bin/python3
"""
This script will generate the input files required for the next steps. 
The input files that will be created are:
1. primer3 file
2. secundary structure prediction file
3. SNP file
4. filter file

The script requires the following arguments:
1. -i: path to the input bed file
2. -p: number of primers
3. -n: difference between primers
4. -a: minimum melting temperature
5. -b: maximum melting temperature
6. -c: optimal melting temperature
7. -d: melting temperature difference
8. -e: minimum GC content
9. -f: maximum GC content
10. -g: optimal GC content
11. -j: minimum amplicon length
12. -k: maximum amplicon length
Note: the maximum amplicon length can be set to 0, in which case it will be calculated based on the sequence length.
"""
####################################################################################################
####################################  1. Parse arguments  ##########################################
####################################################################################################

## import all libraries
import argparse
import math

## get info required to generate the files
parser = argparse.ArgumentParser(description='give arguments to this script for in file generation required for the next steps')
# Input file
parser.add_argument('-i', nargs=1, required=True, help='input fusionRNA bed file, 0-based')
# primer3 parameters
parser.add_argument('-p', nargs=1, required=True, help='nr of primers')
parser.add_argument('-n', nargs=1, required=True, help='nr of difference between primers')
parser.add_argument('-a', nargs=1, required=True, help='min TM')
parser.add_argument('-b', nargs=1, required=True, help='max TM')
parser.add_argument('-c', nargs=1, required=True, help='opt TM')
parser.add_argument('-d', nargs=1, required=True, help='TM diff')
parser.add_argument('-e', nargs=1, required=True, help='min GC')
parser.add_argument('-f', nargs=1, required=True, help='max GC')
parser.add_argument('-g', nargs=1, required=True, help='opt GC')
parser.add_argument('-j', nargs=1, required=True, help='min amp length')
parser.add_argument('-k', nargs=1, required=True, help='max amp length') # this param is set to 0, so that it can be adjusted depending on temp_l if the user does not supply amp_max


## parse the arguments
args = parser.parse_args()
# input file
input_bed = open(args.i[0])
# number of primers
nr = args.p[0]
# difference between primers
diff = args.n[0]

    
####################################################################################################
####################################    2. Processing     ##########################################
####################################################################################################

## read first (and only) line of input bed file
fusionRNA = input_bed.read()

# retrieve chrom start end info
fusion_chrom1 = fusionRNA.split()[0]
fusion_end = int(fusionRNA.split()[1]) # 0-based
fusion_chrom2 = fusionRNA.split()[2]
fusion_start = int(fusionRNA.split()[3]) # 0-based
fusion_ID = fusionRNA.split()[4]


# retrieve sequence
fasta = open("fasta_out.txt")
fasta_track = open("fasta_track.txt")

left = ""
right = ""

# get left and right side of the breakpoint (track file as control)
for seq, seq_type in zip(fasta, fasta_track):
	if seq_type[0:4] == '>end':
		left = seq.rstrip() + left
	elif seq_type[0:6] == '>start':
		right = seq.rstrip() + right

## paste both sides of BP together
sequence = left + right

## get max amp length => by user or depending on temp_length
# temp_length 
length = math.floor(len(sequence) / 2)
amp_max = args.k[0]
if int(amp_max) == 0:
	amp_max = str(int(length + 30))

####################################################################################################
####################################     Primer3 File     ##########################################
####################################################################################################

# sequence based text file (0-based)
output = open("input_primer3_" + fusion_ID + ".txt", "w")
output.write("SEQUENCE_ID=" + fusion_ID + "_" + fusion_chrom1 + "_" + str(fusion_end) + "_*_" + fusion_chrom2 + "_" +str(fusion_start) + "\n")
output.write("SEQUENCE_TEMPLATE=" + sequence + "\n")
output.write("SEQUENCE_TARGET=" + str(int(length-1)) + ",1\nPRIMER_NUM_RETURN=" + nr + "\n")
output.write("PRIMER_MIN_THREE_PRIME_DISTANCE=" + diff + '\n')
output.write("PRIMER_PRODUCT_SIZE_RANGE=" + args.j[0] + '-' + amp_max + '\n')
output.write("PRIMER_MIN_TM=" + args.a[0] + '\n')
output.write("PRIMER_MAX_TM=" + args.b[0]+ '\n')
output.write("PRIMER_OPT_TM=" + args.c[0]+ '\n')
output.write("PRIMER_PAIR_MAX_DIFF_TM=" + args.d[0]+ '\n')
output.write("PRIMER_MIN_GC="  + args.e[0]+ '\n')
output.write("PRIMER_MAX_GC=" + args.f[0]+ '\n')
output.write("PRIMER_OPT_GC_PERCENT="  + args.g[0]+ '\n')
output.close()

####################################################################################################
#######################     secundary structure predicition file     ###############################
####################################################################################################

# fasta file (sequence based)
output = open("input_get_sec_structure_" + fusion_ID + ".fasta", "w")
output.write("> " + fusion_ID + "\n")
output.write(sequence + "\n")
output.close()

####################################################################################################
####################################       SNP file       ##########################################
####################################################################################################

# bed file for SNPs (0-based => use original bed file annotation)
output = open("input_SNP_" + fusion_ID + ".bed", 'w')
output.write(fusion_chrom1+ "\t" + str(fusion_end-length) + "\t" + str(fusion_end) + "\t" + fusion_chrom2 + "\t" + str(fusion_start) + "\t" + str(fusion_start + length) + "\t"+ str(fusion_ID) + "\t" + str(int(length)) + "\n") # 0-based
output.close()

####################################################################################################
####################################     Filter  File     ##########################################
####################################################################################################

# bed file for filtering
output = open("input_filter_" + fusion_ID + ".bed", 'w')
output.write(fusionRNA + '\n')
output.close()