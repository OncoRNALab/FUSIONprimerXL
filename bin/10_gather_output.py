#!/usr/bin/python3
"""
This script uses the filtered primers and annotation file from the sequence retrieval.
It filters the primers according to the filter string and writes the selected primers to a new file.
The first amplicon that passes (which has the shortest amplicon length) is selected.
If no primer pair passes all filters, a message is written to the output file.
If no primers were designed for the fusion transcript, a message is written to the output file.
The requires following input:
	- filtered primers file
	- annotation file
The output is a file with the selected primer.
"""
####################################################################################################
####################################  1. Parse arguments  ##########################################
####################################################################################################

import argparse
import os

parser = argparse.ArgumentParser(description='give arguments to gather output script')
parser.add_argument('-i', nargs=1, required=True, help='file with filtered primers, 0 based')
parser.add_argument('-a', nargs=1, required=True, help='annotation of end and start fusion pos')
args = parser.parse_args()
input_file = args.i[0]

####################################################################################################
#####################################   2. Processing   ############################################
####################################################################################################

# get the primers and fusion ID
all_primers = open(input_file)
fusion_ID = str(input_file).split('_')[3]

# create output file
output = open('selected_primers_' + fusion_ID + '.txt', "w")

# get splice info
splice_info = open(args.a[0])
fusion_id, info = splice_info.readline().rstrip().split('**')

#############################################################################################################

# A. all_primers (filtered_primers) is empty if primer3 was not able to design primers for this fusion transcript
if os.path.getsize(input_file) == 0:
	fusion_ID = str(input_file).split('_')[3]
	chrom1 = str(input_file).split('_')[4]
	end = str(input_file).split('_')[5]
	chrom2 = str(input_file).split('_')[6]
	start = str(input_file).split('_')[7]
	output.write(fusion_ID + '\t' + chrom1 + '\t' + end + '\t' + chrom2 + '\t' + start + "\tprimer3 was not able to design primers for this fusion transcript; try less strict settings\n")

#############################################################################################################

# B. if primer3 was able to design primers for this fusion transcript
else: 
	primer_found = "no"
	# filter primers according to the filter string
	for primer in all_primers:
		filter_str = primer.split()[-1]
		if (primer_found == "no") and (filter_str == "PASS"):
			output.write(primer.rstrip() + "\t" + info.lstrip('/') + "\n")
			primer_found = 'yes'
	# if no primer passed all filters
	if primer_found == "no":
		output.close()
		output = open('selected_primers_' + fusion_ID + '.txt', "w")
		output.write(info.lstrip('/') + "\tno primer pair passed all filters for this fusion transcript\n")

# close files
output.close()
all_primers.close()

		