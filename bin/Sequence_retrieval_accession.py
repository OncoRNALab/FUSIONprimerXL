#!/usr/bin/env python3

#############################################################################################################################################
###############################################   Description of the script   ###############################################################
#############################################################################################################################################
"""
This script imports sequences of specified accession numbers from the NCBI database
Uses the start and end position specified by the user parts of the sequence are joined with * indicating the breakpoint (and of sequence 1,
beginning of sequence 2). The program will use the defined breakpoint and take 150 nucleotides from each side of the breakpoint. If no 150
nucleotides are available the program will take as many nucleotides as possible. These sequences are then written to a file in the output 
directory and contain a header to check if the source is correct

The provided input file should be a csv file with the following format:
email address
Accession 1, Accession2, End of sequence 1, Start of sequence 2

Note: The start and end positions are 1 based counting, meaning that the first position is 1 and the last position is the length of the sequence

Example:
A.N.Other@example.com
NM_000546, NM_000492, 153, 10
NM_001371558, NM_000518, 60, 601

The output file will be a text file with the following format:
> DEFINITION ACCESSION 1: defenition; SOURCE: source
> DEFINITION ACCESSION 2: defenition; SOURCE: source
> Fusion1_NM_000546_NM_000492
sequence1*sequence2

The output file will be saved in the output directory as input_seq.txt
"""
#############################################################################################################################################
######################################################   parameters   #######################################################################
#############################################################################################################################################

# Requires biopython and os
import argparse
import csv
from Bio import Entrez
from Bio import SeqIO
import os

# Argument parser
parser = argparse.ArgumentParser(description='Give a csv file with accessions and positions')
parser.add_argument('-c', '--input_csv', nargs=1, required=True, help='input csv file')
parser.add_argument('-o', '--output_dir', nargs=1, required=True, help='output directory')
args = parser.parse_args()

# save current working directory
Base_dir = os.getcwd()
# check if input file exists
if not os.path.exists(args.input_csv[0]):
    exit("Error: Input file does not exist")
# check if it is a csv file
if not args.input_csv[0].endswith(".csv"):
    exit("Error: Input file is not a csv file")
# check if output directory exists
if not os.path.exists(args.output_dir[0]):
    exit("Error: Output directory does not exist")

#############################################################################################################################################
###############################################   parsing input form the input file   #######################################################
#############################################################################################################################################
# Getting input from the input csv file
line_nr = 0
starts = []
ends = []
Accessions = []

# parse the input from the csv file
with open(args.input_csv[0], 'r') as csvfile:
  reader = csv.reader(csvfile)
  for line in reader:
    # getting the email address
    if line_nr == 0:
        email = str(line[0])
        line_nr += 1
    # getting the accessions and positions
    else:
        line = [element.strip() for element in line]
        Accessions = Accessions + line[0:2]
        ends.append(line[2])
        starts.append(line[3])
csvfile.close()

# NCBI requires an email adress for authentication
# If used excessive an email will be sent to you
Entrez.email = email
del email

# Change the positions from 1 based to 0 based counting (1 in start is from beginning, 0 in end is to end)
start = []
end = []
for i in starts:
    i = int(i)
    if i == 0:
        exit("Error: Start position cannot be 0, please use 1 based counting")
    else:
        i -= 1
    start.append(i)
for i in ends:
    i = int(i)
    if i == 0:
        exit("Error: End position cannot be 0, please use 1 based counting, please use 1 based counting")
    else:
        end.append(i)

# Check if the number of start and end positions are equal
if len(starts) != len(ends):
    exit("Error: The number of start and end positions are not equal " + str(len(starts)) + " start positions and " + str(len(ends)) + " end positions were provided")
del starts
del ends

#############################################################################################################################################
###############################################   Retrieving the sequences   ################################################################
#############################################################################################################################################

# counters 
ID = 0
seq_counter = 0
pos_counter = 0
# move to output directory
os.chdir(args.output_dir[0])

# function to take the surrounding 150 nucleotides of the breakpoint
def left_150(sequence):
    if len(sequence) < 150:
        return sequence
    else:
        return sequence[-150:]
def right_150(sequence):
    if len(sequence) < 150:
        return sequence
    else:
        return sequence[:150]

# Open the sequence output file
sequence_file = open("input_seq.txt", "w")
# Loop through the accessions
for Accession in Accessions:
    Accession = str(Accession)
    try:
        # Fetch the genbank record
        handle = Entrez.efetch(db="nucleotide", id=Accession, rettype="gb", retmode="text")
        sequence_data = handle.read()
        handle.close()
        # write the record to a temporary file
        output_file = open("sequence.gb", "w")
        output_file.write(sequence_data)
        output_file.close()
    except:
        print("Error: Accession could not be found not found")
        print("Tried to access: " + Accession + "\nAt NCBI database nucleotide\nPlease change the accession number and try again")
        exit()
    # parse the genbank record file (nucleotide database)
    input_file = open("sequence.gb", "r")
    sequence_record = SeqIO.read(input_file, "genbank")
        # get the sequence
    sequence = sequence_record.seq
    # First sequence
    if seq_counter % 2 == 0:
        # start and end position of the sequence
        end_pos = end[pos_counter]
        start_pos = start[pos_counter]
        # Slice the sequence
        if end_pos > len(sequence):
            end_pos = len(sequence)+1
        sequence1 = sequence[:end_pos]
        sequence1 = left_150(sequence1)
        # extract the header information
        first_Accession = Accession 
        defenition1 = sequence_record.description
        source1 = sequence_record.annotations.get("source", None)
        # write the first header line
        sequence_file.write("> DEFENITION ACCESSION 1: " + defenition1 + "; SOURCE: " + source1 + "\n") 
    # Second sequence
    else:
        # Slice the sequence
        if start_pos > len(sequence):
            exit("Error: The start position is larger than the sequence length for " + Accession)
        sequence2 = sequence[start_pos:]
        sequence2 = right_150(sequence2)
        # extract the header information
        defenition2 = sequence_record.description
        source2 = sequence_record.annotations.get("source", None)
        # write the second header line
        sequence_file.write("> DEFENITION ACCESSION 2: " + defenition2 + "; SOURCE: " + source2 + "\n")
        # write the third header line
        sequence_file.write("> Fusion" + str() + "_" + str(first_Accession) + "_" + str(Accession) +"\n")
        # write the sequences to the file
        sequence_file.write(str(sequence1))
        sequence_file.write("*")
        sequence_file.write(str(sequence2))
        sequence_file.write("\n")
        # Terminal progress update
        print("Retrieved sequence for: " + first_Accession +" and " + Accession)
        # Increase fusion gene ID
        ID+=1
        pos_counter+=1
    # close the input file (genbank record)
    input_file.close()
    # remove the temporary genbank record file
    os.remove("sequence.gb")
    # increase the sequence counter
    seq_counter+=1
# close the finished output file
sequence_file.close()

# move back to the base directory
os.chdir(Base_dir)
