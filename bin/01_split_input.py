#!/usr/bin/env python3


import argparse

parser = argparse.ArgumentParser(description='give arguments to main primer_xc script')
parser.add_argument('-i', nargs=1, required=True, help='input fusion seq file')
args = parser.parse_args()

input_file = open(args.i[0])
count_lines = 0
for line in input_file:
	count_lines +=1
input_file.close()

fusion_nr = len(str(count_lines))
fusion_nr = "fusion{:0" + str(fusion_nr) +"d}"

# also make a file containing all fusions for the end
all_fusion_file = open('all_fusion.txt', 'w')


input_file = open(args.i[0])

all_fusion = [] # to check if there are any doubles

ID = 0
for fusion in input_file:

	ID_str = fusion_nr.format(ID)

	all_fusion.append(fusion)

	ind_fusion_file = open(ID_str + ".txt", "w")
	ind_fusion_file.write(fusion + '\t' + ID_str + '\n')
	ind_fusion_file.close()

	all_fusion_file.write(ID_str + '\t' + fusion)

	
	ID += 1

all_fusion_file.close()

def checkIfDuplicates(listOfElems):
    ''' Check if given list contains any duplicates '''
    if len(listOfElems) == len(set(listOfElems)):
        return False
    else:
        return True

if checkIfDuplicates(all_fusion):
	raise SystemExit('One or more fusion sequences are present more than once in your input file. Please make sure each fusion sequence is unique.')

input_file.close()