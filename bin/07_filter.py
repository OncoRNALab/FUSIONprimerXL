#!/usr/bin/python3


import argparse
import os

parser = argparse.ArgumentParser(description='give arguments to filter script')
parser.add_argument("-A", nargs = 1, required=True, help='input file with fusion + ID')
parser.add_argument("-P", nargs = 1, required=True, help='intput file with all primers')
parser.add_argument('-t', nargs=1, required=True, help='the output folding file from the template')
parser.add_argument('-a', nargs=1, required=True, help='the output folding file from the amplicon')
parser.add_argument('-b', nargs=1, required=True, help='the output specificity file')
parser.add_argument('-p', nargs=1, required=True, help='specificity filter')

args = parser.parse_args()

all_fusion = open(args.A[0])
spec_filter = args.p[0]

# get general info fusion

fusion = all_fusion.read()

sequence = fusion.split()[0]
fusion_ID = fusion.split()[1]
length = len(sequence)/2


# get specificity info from file
spec = open(args.b[0])
avoid_spec = []
all_lines = []

for line in spec:
	all_lines.append(line)
spec.close()


# rules from https://academic.oup.com/clinchem/article/59/10/1470/5622018?login=true fig 6
for i in range(0, len(all_lines) - 1, 2):
	fwd_spec = all_lines[i].split()
	rev_spec = all_lines[i+1].split()
	if fusion_ID == all_lines[i].split("_")[0]:
		fwd_MM, rev_MM = 0, 0
		if len(fwd_spec) > 7:
			fwd_MM = fwd_spec[7].count('>')
		if len(rev_spec) > 7:
			rev_MM = rev_spec[7].count('>')
		if fwd_MM > 0 or rev_MM > 0:
			if spec_filter == 'strict':
				if fwd_MM + rev_MM < 5:
					avoid_spec.append(all_lines[i].split("_")[2])
			if spec_filter == 'loose':
				if (fwd_MM + rev_MM < 3) or (fwd_MM == 1 & rev_MM == 2) or (fwd_MM == 2 & rev_MM == 1):
					avoid_spec.append(all_lines[i].split("_")[2])
		else:
			avoid_spec.append(all_lines[i].split("_")[2])


# get amplicon folding info into dict
amp_fold = open(args.a[0])
amp_fold_avoid = {}
amp_fold_avoid_pos = {}

for lin in amp_fold:
	primer_ID = str(lin.split()[0].split("_")[1])
	deltag = float(lin.rstrip().split()[1])
	str_pos = lin.rstrip().split('\t')[2]
	amp_fold_avoid[primer_ID] = deltag
	amp_fold_avoid_pos[primer_ID] = str_pos

amp_fold.close()


# filter all primers and write to file

all_primers = open(args.P[0])
primer_file = open("all_primers/filtered_primers_" + fusion_ID + ".txt", "a")

# to make log file
total_primers = 0
passed = 0
failed_spec = 0
#failed_SNP = 0
failed_str_temp = 0
failed_str_amp = 0
design = 1
primer_found = 0

if os.path.getsize(args.P[0]) == 0:
	design = 0


for primer in all_primers:

	total_primers += 1

	primer_ID = str(primer.split()[1])
	FWD_pos = int(primer.split()[4])
	FWD_len = int(primer.split()[5])
	REV_pos = int(primer.split()[6])
	REV_len = int(primer.split()[7])

	FWD_poss = range(FWD_pos, FWD_pos + FWD_len)
	REV_poss = range(REV_pos + 1 - REV_len, REV_pos + 1)

	filter_str = ""

	# check specificity
	if primer_ID in avoid_spec:
		filter_str = filter_str + "FAIL_specificity_"

	# check folding amplicon

	fold_amp_avoid = amp_fold_avoid_pos["primer" + str(primer_ID)]

	if fold_amp_avoid != '[]':
		fold_amp_avoid = fold_amp_avoid.replace('[', "").replace(']', '').split(', ')
		fold_amp_avoid = [ int(x) for x in fold_amp_avoid ]

	print(fold_amp_avoid)

	if amp_fold_avoid["primer" + str(primer_ID)] < -15:
		filter_str = filter_str + 'FAIL_fold_amplicon_'
	elif (amp_fold_avoid["primer" + str(primer_ID)] < -5) & any(x in range(FWD_len + 1) for x in fold_amp_avoid):
		filter_str = filter_str + 'FAIL_fold_amplicon_'
	elif (amp_fold_avoid["primer" + str(primer_ID)] < -5) & any(x in range(REV_pos - REV_len - FWD_pos, REV_pos - FWD_pos + 1) for x in fold_amp_avoid):
		filter_str = filter_str + 'FAIL_fold_amplicon_'



	# if all tests succeded => primer pair passed filters
	if filter_str == "":
		filter_str = "PASS_"
		primer_found = 1

	# gather info to add to log file
	if filter_str == "PASS_":
		passed += 1
	if "FAIL_specificity" in filter_str:
		failed_spec += 1
	if "FAIL_fold_template" in filter_str:
		failed_str_temp += 1
	if "FAIL_fold_amplicon" in filter_str:
		failed_str_amp += 1


	primer_file.write(sequence + "\t" + primer.rstrip() + '\t' + filter_str[0:len(filter_str)-1] + "\n")

primer_file.close()
all_primers.close()

# print log file
log_file = open("log_file_" + fusion_ID + ".txt", "a")
log_file.write(fusion_ID + '\t' + 
	str(design) + "\t" + str(primer_found) + "\t" + str(total_primers) + "\t" + 
	str(passed) + '\t' + str(failed_spec) + "\t" + str(failed_str_amp) + '\n' )
log_file.close()
