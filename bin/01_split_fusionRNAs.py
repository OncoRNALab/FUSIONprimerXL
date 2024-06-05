#!/usr/bin/env python3

import argparse
"""
This script is used to split up the input bed file into different fusion sequences generating a 
bed file for each fusion sequence and a file text file containing all the fusions.
The bed file should be in the following format:
chrom1  end     chrom2  start
example:
chr1	17055   chr1	17232
"""
####################################################################################################
####################################  1. Parse arguments  ##########################################
####################################################################################################
parser = argparse.ArgumentParser(description='Give input to the FUSIONprimerXL pipeline')
# input fusion seq file
parser.add_argument('-i', nargs=1, required=True, help='input fusionRNA bed file, 0-based')
args = parser.parse_args()

input_bed = args.i[0]
####################################################################################################
####################################  2. Processing  ###############################################
####################################################################################################

# open the input bed file and count the number of lines
input_file = open(input_bed)
count_lines = 0
for line in input_file:
	count_lines +=1
input_file.close()

# create a format string for fusionRNA ID
fusion_nr = len(str(count_lines))
fusion_nr = "fusion{:0" + str(fusion_nr) +"d}"

# also make a file containing all fusions for the end
all_fusions_file = open('all_fusions.txt', 'w')

# open the input bed file again
input_file = open(input_bed)

all_fusions = [] # to check if there are any doubles

# Generate the all_fusions file and the individual bed files
ID = 0
for fusionRNA in input_file:

    ID_str = fusion_nr.format(ID)

    chrom1 = fusionRNA.split()[0]
    end= str(fusionRNA.split()[1])
    chrom2 = fusionRNA.split()[2]
    start = str(fusionRNA.split()[3])
    fusion_str = chrom1 + '\t' + end + '\t' + chrom2 + '\t' + start
    # append to the all fusions list
    all_fusions.append(chrom1 + ":" + end + "*" + chrom2 + start)
    # write the fusion to its own bed file
    ind_fusion_file = open(ID_str + ".bed", "w")
    ind_fusion_file.write(fusion_str + '\t' + ID_str + '\n')
    ind_fusion_file.close()
    # write the fusion to the all fusions file
    all_fusions_file.write(ID_str + '\t' + fusion_str + '\n')

    
    ID += 1

all_fusions_file.close()

# check for duplicates
def checkIfDuplicates(listOfElems):
    ''' Check if given list contains any duplicates '''
    if len(listOfElems) == len(set(listOfElems)):
        return False
    else:
        return True

if checkIfDuplicates(all_fusions):
	raise SystemExit('One or more fusionRNAs are present more than once in your input file. Please make sure each fusionRNA is unique.')

input_file.close()