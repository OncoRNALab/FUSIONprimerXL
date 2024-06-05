#!/usr/bin/python3
"""
This script will summarize the run of the pipeline. It will write a summary to the output directory.
The summary will contain the following information:
- the number of fusionRNAs given as input
- the number of fusionRNAs for which primer pairs were designed and passed all filters
- the number of fusionRNAs for which the (spliced) template sequence was too short and no primer pairs could be designed by primer3
- the number of fusionRNAs for which no primer pairs could be designed by primer3
- the number of fusionRNAs for which none of the primer pairs passed all filters
- the number of primer pairs designed and tested
- the number of primer pairs that passed all filters
- the number of primer pairs that failed the specificity filter
- the number of primer pairs that failed the SNP filter
- the number of primer pairs that failed the secondary structure of the template filter
- the number of primer pairs that failed the secondary structure of the amplicon filter
- the start time of the run
- the end time of the run
- the total run time
It requires the following arguments:
- -l: log file
- -s: start time file
- -o: output directory
- -u: upfront filter (yes or no)
- -a: file with all fusionRNAs
"""
####################################################################################################
####################################  1. Parse arguments  ##########################################
####################################################################################################
import pandas as pd
import argparse
from datetime import datetime

parser = argparse.ArgumentParser(description='give arguments to summary script')
parser.add_argument("-l", nargs = 1, required=True, help='log_file')
parser.add_argument('-s', nargs = 1, required=True, help='start time file')
parser.add_argument('-o', nargs = 1, required=True, help='output_dir')
parser.add_argument('-u', nargs = 1, required=True, help='upfront filter')
parser.add_argument('-a', nargs = 1, required=True, help='file with all fusionRNAs')

args = parser.parse_args()

####################################################################################################
####################################   2. Processing   #############################################
####################################################################################################

# log file
log_file = pd.read_csv(open(args.l[0]), sep = "\t", usecols = ['fusion_ID', 'chrom1', 'end',  'chrom2',  'start', 'design', 'primer_found', 'total_primers',  'passed',  'failed_spec', 'failed_SNP', 'failed_str_temp',  'failed_str_amp'])
# start time
start_time = open(args.s[0])
# output directory
out_dir = args.o[0]
# upfront filter
up_f = args.u[0]
# all fusions
all_fusions = open(args.a[0])


summary = open(out_dir + "/summary_run.txt", "a")

#total_nr_fusion = log_file.shape[0]
total_nr_fusion = 0
all_fusions_dict = {}

# count the number of fusions and add to dictionary
for line in all_fusions:
	total_nr_fusion += 1
	# key = fusion ID, value = [chrom1, end, chrom2, start]
	all_fusions_dict[line.rstrip().split('\t')[0]] = line.rstrip().split('\t')[1:]


# calculate the parameters
primer_found = log_file['primer_found'].sum()
fusion_short = total_nr_fusion - log_file.shape[0]
primer_no_design = total_nr_fusion - log_file['design'].sum() - fusion_short
primer_not_found = total_nr_fusion - primer_found - primer_no_design - fusion_short
total_nr_primers = log_file['total_primers'].sum()
primers_passed = log_file['passed'].sum()
failed_spec = log_file['failed_spec'].sum()
failed_snp = log_file['failed_SNP'].sum()
failed_str_temp = log_file['failed_str_temp'].sum()
failed_str_amp = log_file['failed_str_amp'].sum()
# calculate percentages
primer_found_perc = 100 * primer_found/total_nr_fusion
fusion_short_perc = 100 * fusion_short/total_nr_fusion
primer_not_found_perc = 100 * primer_not_found/total_nr_fusion
primer_no_design_perc = 100 * primer_no_design/total_nr_fusion
primers_passed_perc = 100 * primers_passed/ total_nr_primers
failed_spec_perc = 100 * failed_spec / total_nr_primers
failed_snp_perc = 100 * failed_snp / total_nr_primers
failed_str_temp_perc = 100 * failed_str_temp / total_nr_primers
failed_str_amp_perc = 100 * failed_str_amp / total_nr_primers

# write summary
if up_f == "yes":
	summary.write("{total_nr_fusion} fusionRNA(s) were given as input\n\tfor {primer_found} fusionRNA(s) ({primer_found_perc:.1f} %), primer pair(s) were designed and filtered successfully \n\tfor {fusion_short} fusionRNA(s) ({fusion_short_perc:.1f} %), the (spliced) template sequence was too short and no primer pair(s) could be designed by primer3\n\tfor {primer_no_design} fusionRNA(s) ({primer_no_design_perc:.1f} %), no primer pair(s) could be designed by primer3\n\tfor {primer_not_found} fusionRNA(s) ({primer_not_found_perc:.1f} %), none of the primer pair(s) passed all filters\n\nin total\n\t{total_nr_primers} primer pair(s) were designed and tested\n\t{primers_passed} ({primers_passed_perc:.1f} %) primer pairs passed all filters\n\t{failed_spec} ({failed_spec_perc:.1f} %) primer pair(s) failed the specificity filter\n\t{failed_str_amp} ({failed_str_amp_perc:.1f} %) primer pair(s) failed the secondary structure of the amplicon filter\n\nsee log file for more details per fusionRNA\n\n".format(total_nr_fusion = total_nr_fusion, primer_found = primer_found, primer_found_perc = primer_found_perc, primer_no_design = primer_no_design, primer_no_design_perc = primer_no_design_perc, total_nr_primers = total_nr_primers, primers_passed = primers_passed, primers_passed_perc = primers_passed_perc, failed_spec = failed_spec, failed_spec_perc = failed_spec_perc, failed_str_amp = failed_str_amp, failed_str_amp_perc = failed_str_amp_perc, primer_not_found = primer_not_found, primer_not_found_perc = primer_not_found_perc, fusion_short = fusion_short, fusion_short_perc = fusion_short_perc))
else:
	summary.write("{total_nr_fusion} fusionRNA(s) were given as input\n\tfor {primer_found} fusionRNA(s) ({primer_found_perc:.1f} %), primer pair(s) were designed and filtered successfully \n\tfor {fusion_short} fusionRNA(s) ({fusion_short_perc:.1f} %), the (spliced) template sequence was too short and no primer pair(s) could be designed by primer3\n\tfor {primer_no_design} fusionRNA(s) ({primer_no_design_perc:.1f} %), no primer pair(s) could be designed by primer3\n\tfor {primer_not_found_perc} fusionRNA(s) ({primer_not_found:.1f} %), none of the primer pair(s) passed all filters\n\nin total\n\t{total_nr_primers} primer pair(s) were designed and tested\n\t{primers_passed} ({primers_passed_perc:.1f} %) primer pair(s) passed all filters\n\t{failed_spec} ({failed_spec_perc:.1f} %) primer pair(s) failed the specificity filter\n\t{failed_snp} ({failed_snp_perc:.1f} %) primer pair(s) failed the SNP filter\n\t{failed_str_temp} ({failed_str_temp_perc:.1f} %) primer pair(s) failed the secondary structure of the template filter\n\t{failed_str_amp} ({failed_str_amp_perc:.1f} %) primer pair(s) failed the secondary structure of the amplicon filter\n\nsee log file for more details per fusionRNA\n\n".format(total_nr_fusion = total_nr_fusion, primer_found = primer_found, primer_found_perc = primer_found_perc, primer_no_design = primer_no_design, primer_no_design_perc = primer_no_design_perc, total_nr_primers = total_nr_primers, primers_passed = primers_passed, primers_passed_perc = primers_passed_perc, failed_spec = failed_spec, failed_spec_perc = failed_spec_perc, failed_snp = failed_snp, failed_snp_perc = failed_snp_perc, failed_str_temp = failed_str_temp, failed_str_temp_perc = failed_str_temp_perc, failed_str_amp = failed_str_amp, failed_str_amp_perc = failed_str_amp_perc, primer_not_found = primer_not_found, primer_not_found_perc = primer_not_found_perc, fusion_short = fusion_short, fusion_short_perc = fusion_short_perc))

# write start and end time
start_str = start_time.read().rstrip()
start = datetime.strptime(start_str, "%d/%m/%Y %H:%M:%S")
end = datetime.now()
end_str = end.strftime("%d/%m/%Y %H:%M:%S")

diff = int((end - start).total_seconds())

summary.write("run started at " + start_str + " (UTC)\nrun ended at " + end_str + " (UTC)\n" + "total run time: " + str(diff) + " seconds\n")

summary.close()


# also add fusionRNAs that are not in log and filtered primers
log_file = open(args.l[0])

log_file_list = []
for line in log_file:
	log_file_list.append(line.split('\t')[0])

log_file.close()
log_file = open(args.l[0], 'a')
filtered = open('filtered_primers.txt', 'a')

for ID in all_fusions_dict:
	if ID not in log_file_list:
		log_file.write(ID + ('\tthis (spliced) fusionRNA template sequence is too short to design primers for\n'))
		filtered.write(ID + ('\tthis (spliced) fusionRNA template sequence is too short to design primers for\n'))

log_file.close()
filtered.close()
