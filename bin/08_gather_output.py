#!/usr/bin/python3

import argparse
import os

parser = argparse.ArgumentParser(description='give arguments to main primer_xc script')
parser.add_argument('-i', nargs=1, required=True, help='input fusion seq file')
args = parser.parse_args()
input_file = args.i[0]

all_primers = open(input_file)
fusion_ID = str(input_file).split('_')[3]
output = open('selected_primers_' + fusion_ID + '.txt', "w")

# get splice info

if os.path.getsize(input_file) == 0:
	fusion_ID = str(input_file).split('_')[3]

	output.write(fusion_ID + "\tprimer3 was not able to design primers for this fusion sequence; try less strict settings\n")

else: 
	primer_found = "no"
	# filter primers according to the filter string
	for primer in all_primers:
		filter_str = primer.split()[-1]
		if (primer_found == "no") and (filter_str == "PASS"):
			output.write(primer.rstrip() + "\t" + "\n")
			primer_found = 'yes'
	if primer_found == "no":
		output.write(('\t'.join(primer.split()[0:2])) + "\tno primer pair passed all filters for this fusion sequence\n")


output.close()

		