#!/usr/bin/python3
"""
This script filters primers based on specificity, SNPs, folding template and folding amplicon.
It takes the following arguments:
-A: input file with fusionRNA + ID (input_filter_fusionID.bed)
-P: intput file with all primers
-l: the nr of nucleotides surrounding the BP at each side
-s: the output snp file
-t: the output folding file from the template
-a: the output folding file from the amplicon
-b: the output specificity file
-p: specificity filter (strict, loose)
-f: filtering SNPs, should be strict or loose (off, strict, loose)
And outputs a file with filtered primers and a log file with the results of the filtering proces
The log file contains the following columns:
fusion_ID, chrom1, end, chrom2, start, design, primer_found, total_primers, passed, failed_spec, failed_SNP, failed_str_temp, failed_str_amp
"""
####################################################################################################
####################################  1. Parse arguments  ##########################################
####################################################################################################
import argparse
import os

parser = argparse.ArgumentParser(description='give arguments to filter script')
parser.add_argument("-A", nargs = 1, required=True, help='input file with fusionRNA + ID (input_filter_fusionID.bed)')
parser.add_argument("-P", nargs = 1, required=True, help='intput file with all primers')
parser.add_argument('-l', nargs=1, required=True, help='the nr of nucleotides surrounding the BP at each side')
parser.add_argument('-s', nargs=1, required=True, help='the output snp file')
parser.add_argument('-t', nargs=1, required=True, help='the output folding file from the template')
parser.add_argument('-a', nargs=1, required=True, help='the output folding file from the amplicon')
parser.add_argument('-b', nargs=1, required=True, help='the output specificity file')
parser.add_argument('-p', nargs=1,choices=['strict','loose'], required=True, help='specificity filter (strict, loose)')
parser.add_argument('-f', nargs=1,choices=['off','strict','loose'], required=True, help='filtering SNPs, should be strict or loose (off, strict, loose)')

args = parser.parse_args()

all_fusion = open(args.A[0])
length = int(args.l[0])
SNP_filter = args.f[0]
spec_filter = args.p[0]

####################################################################################################
#####################################    2. Processing   ###########################################
####################################################################################################

# get general info fusionRNA
fusionRNA = all_fusion.read()

chrom1 = fusionRNA.split()[0]
end = int(fusionRNA.split()[1])
chrom2 = fusionRNA.split()[2]
start = int(fusionRNA.split()[3])
fusion_ID = fusionRNA.split()[4]



# get specificity info from file (bowtie output)
spec = open(args.b[0])
avoid_spec = []
all_lines = []

for line in spec:
	all_lines.append(line)
spec.close()

####################################################################################################
#################################    A. specificity filter  ########################################
####################################################################################################

# rules from https://academic.oup.com/clinchem/article/59/10/1470/5622018?login=true fig 6 are used as filter criteria
# The criteria are:
# Strict: 
# - no mismatches in either primer = off-target and is discarted
# - primers with at least 4 mismatches for a single primer = no off-target and is kept
# - primers with a total of at least 5 mismatches between both primers = no off-target and is kept
# Loose:
# - primers with at least 3 mismatches for a single primer = no off-target and is kept
# - primers with a total of at least 4 mismatches between both primers = no off-target and is kept

# MM primer1	MM primer2	sum MM	loose			strict
# 0				0			0		off-target		off-target
# 1				0			1		off-target		off-target
# 0				1			1		off-target		off-target
# 2				0			2		off-target		off-target
# 0				2			2		off-target		off-target
# 1				1			2		off-target		off-target
# 2				1			3		off-target		off-target
# 1				2			3		off-target		off-target
# 3				0			3		no off-target	off-target
# 0				3			3		no off-target	off-target
# 1				3			4		no off-target	off-target
# 3				1			4		no off-target	off-target
# 2				2			4		no off-target	off-target
# 4				0			4		no off-target	no off-target
# 0				4			4		no off-target	no off-target
# 2				3			5		no off-target	no off-target
# 3				2			5		no off-target	no off-target
# 4				1			5		no off-target	no off-target
# 1				4			5		no off-target	no off-target
# 3				3			6		no off-target	no off-target

# get the specificity info from the file
for i in range(0, len(all_lines) - 1, 2):
	fwd_spec = all_lines[i].split()
	rev_spec = all_lines[i+1].split()
	# loop
	if fusion_ID == all_lines[i].split("_")[0]:
		# get the number of mismatches
		fwd_MM, rev_MM = 0, 0
		if len(fwd_spec) > 7:
			fwd_MM = fwd_spec[7].count('>')

		if len(rev_spec) > 7:
			rev_MM = rev_spec[7].count('>')

		# check the number of mismatches
		if fwd_MM > 0 or rev_MM > 0:
			# STRICT
			if spec_filter == 'strict':
				# at least 5 combined mismatches
				if (fwd_MM + rev_MM < 5):
					# at least 4 mismatches in one primer
					if (fwd_MM == 4 and rev_MM == 0) or (fwd_MM == 0 and rev_MM == 4) :
						pass
					# reject
					else:
						avoid_spec.append(all_lines[i].split("_")[2])
			# LOOSE
			if spec_filter == 'loose':
        		# at least 3 mismatches in one primer
				if (fwd_MM < 3 and rev_MM < 3):
					# at least 4 combined mismatches
					if (fwd_MM == 2 and rev_MM == 2):
						pass
					# reject
					else:
						avoid_spec.append(all_lines[i].split("_")[2])
		# no mismatches: off-target
		else:
			avoid_spec.append(all_lines[i].split("_")[2])


####################################################################################################
#################################    3. masked positions    ########################################
####################################################################################################

# A. Template secundary structure
## get folding template from file
fold_temp = open(args.t[0])

skip = fold_temp.readline()
skip = fold_temp.readline()
skip = fold_temp.readline()
fold_temp_avoid = fold_temp.readline()

## make a list from the positions to avoid based on folding template
if fold_temp_avoid != '[]':
	fold_temp_avoid = fold_temp_avoid.replace('[', "").replace(']', '').split(', ')
	fold_temp_avoid = [ int(x) for x in fold_temp_avoid ]


# B. Single nucleotide polymorphisms
## get SNP positions from file
SNPs = open(args.s[0]).readline()
## if no SNPs are present, set SNP_avoid to none
SNP_avoid = 'none'

## make a list from the positions
if SNPs != '[]':
	SNPs = SNPs.replace('[', "").replace(']', '').split(', ')
	SNP_avoid = [ int(x) for x in SNPs ]



# C. Amplicon secundary structure
## get amplicon folding info into dicts
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

####################################################################################################
#################################    4. LOG file generation    #####################################
####################################################################################################

# filter all primers and write to file
# a filter string will be used to check if the primer passes certain conditions or not

# create an output file in the all primers folder
all_primers = open(args.P[0])
primer_file = open("all_primers/all_primers_" + fusion_ID + '_' + chrom1 + '_' + str(end) + '_' + chrom2 + '_' + str(start) + "_.txt", "a")

# to make log file
total_primers = 0
passed = 0
failed_spec = 0
failed_SNP = 0
failed_str_temp = 0
failed_str_amp = 0
design = 1
primer_found = 0

# check if file is empty
if os.path.getsize(args.P[0]) == 0:
	design = 0

# loop over all primers
for primer in all_primers:
	# add 1 to primer count
	total_primers += 1
	# get info from primer
	primer_ID = str(primer.split()[1])
	FWD_pos = int(primer.split()[4])
	FWD_len = int(primer.split()[5])
	REV_pos = int(primer.split()[6])
	REV_len = int(primer.split()[7])

	# create ranges from the forward and reverse positions
	FWD_poss = range(FWD_pos, FWD_pos + FWD_len)
	REV_poss = range(REV_pos - REV_len, REV_pos)
	# create filter string
	filter_str = ""

	# check specificity (see earlier)
	if primer_ID in avoid_spec:
		filter_str = filter_str + "FAIL_specificity_"

	####################################################################################################
	#########################################   B. SNP filter   ########################################
	####################################################################################################
	if (SNP_filter != 'off') & (SNP_avoid != 'none'):
		# STRICT
		if SNP_filter == 'strict':
			if any(x in FWD_poss for x in SNP_avoid):
				filter_str = filter_str + "FAIL_snp_F_"
			if any(x in REV_poss for x in SNP_avoid):
				filter_str = filter_str + 'FAIL_snp_R_'
		
		# LOOSE
		if SNP_filter == 'loose':
			# 5 end based
			FWD_poss_5end = FWD_poss[0:int(len(FWD_poss)/2)]
			REV_poss_5end = REV_poss[int(len(REV_poss)/2):]

			if any(x in FWD_poss_5end for x in SNP_avoid):
				filter_str = filter_str + "FAIL_snp_F_"

			if any(x in REV_poss_5end for x in SNP_avoid):
				filter_str = filter_str + 'FAIL_snp_R_'

	####################################################################################################
	#################################    C. template sec struc filter  #################################
	####################################################################################################
	# check folding template
	if any(x in FWD_poss for x in fold_temp_avoid):
		filter_str = filter_str + 'FAIL_fold_template_F_'
	if any(x in REV_poss for x in fold_temp_avoid):
		filter_str = filter_str + 'FAIL_fold_template_R_'


	####################################################################################################
	#################################    D. amplicon sec struc filter  #################################
	####################################################################################################
	fold_amp_avoid = amp_fold_avoid_pos["primer" + str(primer_ID)]

	if fold_amp_avoid != '[]':
		# generate list
		fold_amp_avoid = fold_amp_avoid.replace('[', "").replace(']', '').split(', ')
		fold_amp_avoid = [ int(x) for x in fold_amp_avoid ]

	# deltaG check
	## less than -15 is not allowed
	## less than -5 is not allowed if the position of a secundary structure is in the range of the primer
	## greater than -5 is allowed
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
	if "FAIL_snp" in filter_str:
		failed_SNP += 1
	if "FAIL_fold_template" in filter_str:
		failed_str_temp += 1
	if "FAIL_fold_amplicon" in filter_str:
		failed_str_amp += 1

	primer_file.write(primer.rstrip() + "\t" + filter_str[0:len(filter_str)-1] + "\n")

# end of primer loop

primer_file.close()
all_primers.close()

# write log file
log_file = open("log_file_" + fusion_ID + ".txt", "a")
log_file.write(fusion_ID + '\t' + chrom1 + '\t' + str(end) + '\t' + chrom2 + '\t' + str(start) + '\t' + 
	str(design) + "\t" + str(primer_found) + "\t" + str(total_primers) + "\t" + 
	str(passed) + '\t' + str(failed_spec) + "\t" + str(failed_SNP) + "\t" + 
	str(failed_str_temp) + "\t" + str(failed_str_amp) + '\n' )
log_file.close()
