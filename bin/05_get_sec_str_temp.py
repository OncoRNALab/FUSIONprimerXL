#!/usr/bin/python3
"""
This script takes a fasta file as input and returns the minimal free energy structure and delta G in the output file.
From this structure, it also returns the list of indexes that are not allowed to be used as primers because they are in the secondary structure.
This script is used to get the secondary structure of the fusion gene at 60 degrees.
With salt concentration of 0.05 M and DNA nucleic acid type.
The required input:
	-i: input file with the sequence of the fusion gene in fasta format

The script will also create a fusion_ID_ss.ps file that can be used to visualize the secondary structure of the fusion gene.
Example: $gv fusion0_ss.ps to visualize the secondary structure of the fusion gene.
These files are not stored in the output by default when using the pipeline.
"""
####################################################################################################
####################################  1. Parse arguments  ##########################################
####################################################################################################
import os
import argparse

parser = argparse.ArgumentParser(description='give arguments to main primer_xc script')
parser.add_argument('-i', nargs=1, required=True, help='get sec structure input file (fasta file)')
args = parser.parse_args()
input_file = args.i[0]

####################################################################################################
####################################     2. Minimal free energy   ##################################
####################################################################################################

# read in sequence
sequence_file = open(input_file)
fusion_ID = sequence_file.readline().rstrip().replace('> ', '')
sequence = sequence_file.readline().rstrip()

# get the minimal free energy structure
## change temperature to 60 degrees
## change salt to 0.05 M
## change the type of nucleic acid to DNA
command = "RNAfold -p -i " + input_file + " --noconv -T 60.0 -P DNA --salt 0.05"
process = os.popen(command)
output = process.read()
process.close()  # Ensure proper resource management
# Extract stdout from the command
arguments = output.split("\n")

# create an output file
output = open('output_RNAfold_temp_' + fusion_ID + '.txt', 'w')
## fusion_ID
output.write(arguments[0] + "\n")
## structure
structure = arguments[4].split(" ")[0]
output.write(structure + "\n")
## delta_g
delta_g = arguments[3].split(" ",1)[1]
output.write(delta_g+ "\n") # delta_g


####################################################################################################
####################################     3. positions to avoid    ##################################
####################################################################################################

# make list of sec str to avoid
str_not_ok = []
index = 0

for sign in structure:
	if sign == '(' or sign == ')':
		str_not_ok.append(index)
	index += 1


output.write(str(str_not_ok))
output.close()