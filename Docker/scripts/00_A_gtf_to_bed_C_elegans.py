#!/usr/bin/python3

"""
WARNING: this works for gtf file downloaded from ensembl, see README.md, it does not work well for gtf files downloaded from Gencode or the UCSC table browser.

This script will transform a gtf file containing all known genes and their exons into a bed file containing all exons. The bed file will contain the following columns:
1. chromosome
2. start position
3. end position
4. transcript_id_transcript_version_exon_number

The script requires two arguments:
1. -i: path to the input gtf file
2. -o: path to the output bed file
"""
####################################################################################################
####################################  1. Parse arguments  ##########################################
####################################################################################################

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', nargs=1, required=True, help='input gtf file')
parser.add_argument('-o', nargs=1, required=True, help='output bed file')

args = parser.parse_args()

####################################################################################################
###################################  2. process the files  #########################################
####################################################################################################

# open the input gtf file and the output bed file
gtf = open(args.i[0])
bed = open(args.o[0], 'w')

# parse the gtf file and write the bed file
for line in gtf:
	# skip header
	if line[0] != '#': # skip header
		# split the line into columns
		chrom, WormBase, annotation, start, end, point, strand, point2, info = line.split('\t')
		
		# only keep the lines that contain exons
		if annotation == "exon":
			gene_id, transcript_id, exon_number, gene_name, skip= info.split(';',4)

			# extract the relevant information
			transcript_id = transcript_id.lstrip().rstrip().replace('transcript_id ', '').replace('"', '')
			exon_number = exon_number.lstrip().rstrip().replace('exon_number ', '').replace('"', '')
			# write the bed file
			## if the gtf file already contains the chr 
			if chrom[0:3] == "chr": 
				bed.write(chrom + '\t' + str(int(start) - 1 ) + '\t' + end + '\t' + transcript_id + '_exon_' + exon_number + '\n')
			## if the gtf file does not yet contain the chr
			else: 
				bed.write("chr" + chrom + '\t' + str(int(start) - 1 ) + '\t' + end + '\t' + transcript_id + '_' + '_exon_' + exon_number + '\n')

# close the files
gtf.close()
bed.close()
