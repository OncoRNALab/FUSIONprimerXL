#!/usr/bin/python3
"""
This script will validate the BED file that was generated with the script 00_A_gtf_to_bed.py
The script requires two arguments:
1. -i: path to the input bed file
2. -c: path to the chrom.sizes file

The script will check if the chromosome names in the bed file correspond to the chromosome names in the chrom.sizes file.
It will also check if the start and end coordinates are within the chromosome size.
It will output the lines that conain incorrect chromosome names.
"""
####################################################################################################
####################################  1. Parse arguments  ##########################################
####################################################################################################
import argparse


parser = argparse.ArgumentParser(description='give arguments to validate bed script')
parser.add_argument('-i', nargs=1, required=True, help='input fusionRNA bed file, 0-based')
parser.add_argument('-c', nargs=1, required=True, help='chrom file')
args = parser.parse_args()

####################################################################################################
###################################  2. process the files  #########################################
####################################################################################################

# Open the input bed file and the chrom.sizes file
input_bed = open(args.i[0])
chrom_file = open(args.c[0])

# Parse the chrom.sizes file
chrom_sizes = {}

try:
	for line in chrom_file:
		c = line.rstrip().split("\t")
		chrom_sizes[c[0]] = int(c[1])
except:
	print("Error parsing the chrom.sizes file")
	raise

# Validate the bed file
line_number = 0

error = ""


for line in input_bed:
	line_number+=1
	l = line.rstrip().split("\t")

	### check chrom
	chrom = l[0]
	try:
		chromosome_size = chrom_sizes[chrom]
	except:
		error = error + ("Error in bed input file: faulty chromosome name at line %d \n" % (line_number))
		continue

	### check start
	start = l[1]
	try:
		start = int(start)
	except:
		error = error + ("Error in bed input file: parsing start coordinate at line %d \n" % (line_number))
		continue

	if start >= chromosome_size:
		error = error + ("Error in bed input file: start coordinate %d larger than chromosome size %d at line %d \n" % (start, chromosome_size, line_number))
		continue


	### check end
	end = l[2]
	try:
		end = int(end)
	except:
		error = error + ("Error in bed input file: parsing end coordinate at line %d \n" % (line_number))
		continue
	
	if end > chromosome_size:
		error = error + ("Error in bed input file: end coordinate %d larger than chromosome size %d at line %d \n" % (end, chromosome_size, line_number))
		continue

	### check start smaller than start
	if start >= end:
		error = error + ("Error in bed input file: start coordinate %d larger than (or equal to) end coordinate %d at line %d \n" % (start, end, line_number))
		continue

if error != "":
	raise SystemExit(error)

# close the files
input_bed.close()
chrom_file.close()
