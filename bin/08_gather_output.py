#!/usr/bin/python3

import argparse
import os

parser = argparse.ArgumentParser(description='give arguments to main primer_xc script')
parser.add_argument('-i', nargs=1, required=True, help='input circRNA bed file, 0-based')
args = parser.parse_args()
input_file = args.i[0]

all_primers = open(input_file)
circ_ID = str(input_file).split('_')[3]
output = open('selected_primers_' + circ_ID + '.txt', "w")

# get splice info

if os.path.getsize(input_file) == 0:
	circ_ID = str(input_file).split('_')[3]

	output.write(circ_ID + "\tprimer3 was not able to design primers for this circRNA; try less strict settings\n")

else: 
	primer_found = "no"
	# filter primers according to the filter string
	for primer in all_primers:
		filter_str = primer.split()[-1]
		if (primer_found == "no") and (filter_str == "PASS"):
			output.write(primer.rstrip() + "\t" + "\n")
			primer_found = 'yes'
	if primer_found == "no":
		output.write(('\t'.join(primer.split()[0:2])) + "\tno primer pair passed all filters for this circRNA\n")


output.close()

		