#!/usr/bin/python3
"""
This script will add a new line to the input_primer3 file with the positions that are not allowed to be used as primers.
Example of the added line:
SEQUENCE_EXCLUDED_REGION=5,4 10,1 12,6 30,3 34,7 145,6 166,6 179,3 183,7 203,7 211,4 42,1
each pair of numbers represents the position and the length of the region that is not allowed to be used as a primer.
The required input:
	-i: input file with the generated input_primer3 file
	-a: input file with the folding information
	-s: input file with the SNPs
	-f: filter : no/snp/str/yes(both snp and str)
	-l: length of the region surrounding the BP
output:
	- primer3_file_fusion_ID.txt
"""
####################################################################################################
####################################  1. Parse arguments  ##########################################
####################################################################################################

import argparse
# arguments
parser = argparse.ArgumentParser(description='give arguments to filter script')
parser.add_argument("-i", nargs = 1, required=True, help='original primer3 file')
parser.add_argument("-a", nargs = 1, required=True, help='template folding to avoid')
parser.add_argument("-s", nargs=1, required=True, help='snp out')
parser.add_argument("-f", nargs=1, required=True, help='shoud this step be included? no/yes/snp/str')
parser.add_argument('-l', nargs=1, required=True, help='the nr of nucleotides surrounding the BP at each side')

args = parser.parse_args()

####################################################################################################
################################  2. copy the original input file  #################################
####################################################################################################

# length of the region surrounding the BP
length = int(args.l[0])
# original primer3 input file
primer3_file = open(args.i[0])
# get sequence ID info
fusion_ID, chrom1, end,symbol,chrom2,start = primer3_file.readline().replace("SEQUENCE_ID=", "").split("_") # fusion0_chr1_17055_*_chr1_17525
# close and reopen to read from the beginning
primer3_file.close()

# filter on or off
filter_on = args.f[0]
primer3_file = open(args.i[0])
# new primer3 file
primer3_file_new = open("primer3_file_" + fusion_ID + ".txt", 'w')
# copy complete file
for line in primer3_file:
	primer3_file_new.write(line)
	# save length info
	if line[0:17] == "SEQUENCE_TEMPLATE":
		length = len(line.split("=")[1]) / 2

####################################################################################################
####################################  3. structure filter  #########################################
####################################################################################################
avoid_range = []
# if this filter is on (yes meaning both str and snp)
if filter_on == 'yes' or filter_on == "str":
	folding_file = open(args.a[0])

	# get folding info
	skip = folding_file.readline() # ID
	skip = folding_file.readline() # dot bracket notation
	skip = folding_file.readline() # 
	fold_temp_avoid = folding_file.readline()
	folding_file.close()
	

	# put folding info into format for primer 3 (space-separated: position,length)
	## only if there are positions to avoid
	if fold_temp_avoid != '[]':
		fold_temp_avoid = fold_temp_avoid.replace('[', "").replace(']', '').split(', ')
		fold_temp_avoid = [int(x) for x in fold_temp_avoid]

		begin_region = fold_temp_avoid[0]
		index = 0
		length = 1

		# go over all positions and add them to avoid_range
		for position in fold_temp_avoid:
			if index + 1 < len(fold_temp_avoid):
				if position == fold_temp_avoid[index + 1] - 1:
					length += 1
					index += 1
				else:
					avoid_range.append([begin_region, length])
					begin_region = fold_temp_avoid[index + 1]
					index += 1
					length = 1
			else:
				avoid_range.append([begin_region, length])
				length = 1

####################################################################################################
####################################   4. SNP filter   #############################################
####################################################################################################
if filter_on == 'yes' or filter_on == "snp":

	# get SNP info
	SNPs = open(args.s[0]).readline()

	# if there are SNPs to avoid add them to avoid_range
	if SNPs != '[]':
		SNPs = SNPs.replace('[', "").replace(']', '').split(', ')
		SNPs = [ int(x) for x in SNPs ]

		# range is always 1 no need for extending the range
		for snp in SNPs:
			avoid_range.append([snp, 1])

# if filter is on add the avoid_range to the new primer3 file
## else original file is copied with last = line added
if filter_on != "no":
		avoid_range_str = "SEQUENCE_EXCLUDED_REGION="

		for element in avoid_range:
			avoid_range_str = avoid_range_str + str(element[0]) + "," + str(element[1]) + ' '

		primer3_file_new.write(avoid_range_str.rstrip() + '\n')


# add last line with '=' again
primer3_file_new.write("=\n")

# close all files
primer3_file.close()
primer3_file_new.close()
