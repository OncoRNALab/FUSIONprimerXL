#!/usr/bin/env python3


import argparse

parser = argparse.ArgumentParser(description='give arguments to main primer_xc script')
parser.add_argument('-i', nargs=1, required=True, help='input circRNA bed file, 0-based')
args = parser.parse_args()
input_bed = args.i[0]

input_file = open(input_bed)
count_lines = 0
for line in input_file:
	count_lines +=1
input_file.close()

circ_nr = len(str(count_lines))
circ_nr = "circ{:0" + str(circ_nr) +"d}"

# also make a file containing all circs for the end
all_circ_file = open('all_circ.txt', 'w')


input_file = open(input_bed)

all_circ = [] # to check if there are any doubles

ID = 0
for circRNA in input_file:

	ID_str = circ_nr.format(ID)

	all_circ.append(circRNA)

	ind_circ_file = open(ID_str + ".bed", "w")
	ind_circ_file.write(circRNA + '\t' + ID_str + '\n')
	ind_circ_file.close()

	all_circ_file.write(ID_str + '\t' + circRNA)

	
	ID += 1

all_circ_file.close()

def checkIfDuplicates(listOfElems):
    ''' Check if given list contains any duplicates '''
    if len(listOfElems) == len(set(listOfElems)):
        return False
    else:
        return True

if checkIfDuplicates(all_circ):
	raise SystemExit('One or more circRNAs are present more than once in your input file. Please make sure each circRNA is unique.')

input_file.close()