#!/usr/bin/python3
"""
This script will allow the pipeline to run primer3 on the input file and generate the output file with the correct name.
The script will take the input file and primer3 settings file as arguments and run primer3 with the correct output file name.
The output file will be named as output_primer3_fusionID.txt where fusionID is the fusion ID extracted from the input file name.

Input:
    -i: input primer3_file_fusionID.txt
    -x: primer3 settings file
Output:
    output_primer3_fusionID.txt
"""
####################################################################################################
####################################  1. Parse arguments  ##########################################
####################################################################################################
import os
import argparse

parser = argparse.ArgumentParser(description='give the arguments to the primer3 script to execute primer3 with correct output file names')
parser.add_argument('-i',nargs=1, required=True, help='input primer3_file_fusionID.txt')
parser.add_argument('-x', nargs=1, required=True, help='primer3 settings file')
args = parser.parse_args()

input_file = args.i[0]
primer3_settings = args.x[0]

####################################################################################################
####################################     2. Run primer3   ##########################################
####################################################################################################

# Extract the fusion ID from the input file name
fusion_ID = input_file.split('_')[2].split('.')[0]

# construct the command to run primer3
command = "primer3_core --output=output_primer3_" + fusion_ID + ".txt --p3_settings_file=" + primer3_settings + " " + input_file
print(command)

# Run primer3
try:
    os.system(command)
except:
    exit("Error: primer3 failed to run.")




