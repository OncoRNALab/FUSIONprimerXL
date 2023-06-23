#!/usr/bin/python3

from nupack import *
import argparse

parser = argparse.ArgumentParser(description='give arguments to main primer_xc script')
parser.add_argument('-i', nargs=1, required=True, help='nupack input file')
parser.add_argument('-f', nargs=1, required=True, help='is filter on')
args = parser.parse_args()
input_file = args.i[0]


# read in sequence
sequence_file = open(input_file)
fusion_ID = sequence_file.readline().rstrip().replace('> ', '')
sequence = []
sequence.append(sequence_file.readline().rstrip())

# make ouput file
output = open('output_NUPACK_temp_' + fusion_ID + '.txt', 'a')
output.write(fusion_ID + '\n')

if args.f[0] == "on":
	# make model
	my_model = Model(material='DNA', celsius=60, sodium=0.05, magnesium=0.003)

	# MFE proxy structure and energy
	my_mfe = mfe(strands=sequence, model=my_model)
	delta_g = float(my_mfe[0].energy)
	structure = str(my_mfe[0].structure)

	# Print out components of the result for the given complex


	output.write(structure)
	output.write('\n' + str(delta_g) + '\n')

	# make list of sec str to avoid

	str_not_ok = []
	index = 0

	for sign in structure:
		if sign == '(' or sign == ')':
			str_not_ok.append(index)
		index += 1
	output.write(str(str_not_ok))


output.close()
