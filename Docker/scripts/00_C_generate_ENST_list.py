#!/usr/bin/python3

"""
This script will generate a list of ENST identifiers from a MANE file.
The script requires two arguments:
1. -i: path to the input MANE file
2. -c: path to the output ENST list file (name: ENST_list_GRCh38.txt)
"""

####################################################################################################
####################################  1. Parse arguments  ##########################################
####################################################################################################
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', nargs=1, required=True, help='input MANE file')
parser.add_argument('-o', nargs=1, required=True, help='output ENST list file')

args = parser.parse_args() 

####################################################################################################
###################################  2. process the file  ##########################################
####################################################################################################
# open the files
mane = open(args.i[0])
out = open(args.o[0], 'w')

# parse the MANE file
out_list = []

for line in mane:
	if line[0] != '#':
		chrom, BestRefSeq, tr_type, start, end, point, strand, point2, info = line.split('\t')

		if tr_type == "exon":
			gene_id, transcript_id, skip = info.split(';', 2)
			out_list.append(transcript_id.lstrip().rstrip().replace('transcript_id', '').replace('"', ''))
		

# Maintain only unique elements
out_list = set(out_list)

# write the output file
for element in out_list:
	out.write(element[1:] + '\n')

# close the files
mane.close()
out.close()