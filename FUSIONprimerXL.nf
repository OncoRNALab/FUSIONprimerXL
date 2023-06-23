#!/usr/bin/env nextflow

// set default parameters

params.index_bowtie = "$baseDir/assets/GRCh38/index_bowtie"
params.index_bowtie_name = "GRCh38_dna"
params.primer_settings = "$baseDir/assets/primer3plus_settings.txt"
params.input_bed = "example/path"

params.primer3_diff = 1
params.primer3_nr = 20
params.min_tm = 58
params.max_tm = 60
params.opt_tm = 59
params.diff_tm = 2
params.min_gc = 30
params.max_gc = 80
params.opt_gc = 50
params.amp_min = 50
params.amp_max = 0 // this param is set to 0, so that it can be adjusted depending on temp_l if the user does not supply amp_max
params.snp_filter = 'strict'
params.snp_url = 'http://hgdownload.soe.ucsc.edu/gbdb/hg38/snp/dbSnp153Common.bb'
params.temp_str_filter = 'on'
params.spec_filter = 'strict'

params.help = false

// required parameters
input_bed = file(params.input_bed)
index_bowtie = file(params.index_bowtie)

// help message

def helpMessage() {
	log.info"""
	
	Usage:
	
	The typical command for running the pipeline is as follows:
	nextflow run CIRCprimerXL.nf -profile singularity
	
	Mandatory nextflow arguments:
	-profile 		set to 'local' when running locally, set to 'singularity' when running on the HPC

	Mandatory pipeline arguments:
	--input_bed		path to input file with circRNAs in bed format (0-based annotation)
	--index_bowtie		path to bowtie genome index directory
	--index_bowtie_name	the basename of the Bowtie index to be searched (the name of any of the index files up to but not including the final .1.ebwt / .rev.1.ebwt / ...)


	Optional pipeline arguments:
	--primer_settings	path to file with primer3plus settings (see primer3 manual)
	--primer3_diff		the minimum number of base pairs between the 3' ends of any two left primers (see also primer3 PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE)
	--primer3_nr		the number of primers designed by primer3; caution: setting this parameter to a large value will increase running time
	--min_tm		minimum melt temperature of the primers (default: 58)
	--max_tm		maximum melt temperature of the primers(default: 60)
	--opt_tm		optimal melt temperature of the primers(default: 59)
	--diff_tm		maximum difference in melt temperature between the primers(default: 2)
	--min_gc		minimum GC contect of the primers (default: 30)
	--max_gc		maximum GC contect of the primers(default: 80)
	--opt_gc		optimal GC contect of the primers(default: 50)
	--amp_min		minimum amplicon length (default: 60)
	--amp_max		maximum amplicon length (default: 0)
	--snp_filter		when set to 'strict', no common SNPs are allowed in primer sequence; when set to 'loose', common SNPs are allowed in 5' first half of primer
	--snp_url		when using a differente species than human, the correct SNP database url should be provided; alternatively, this paramater can be set to 'off' if no SNP database is available
	--temp_str_filter
	--spec_filter		when set to 'strict', only 2MM + 3MM are allowed; when set to 'loose', 2MM + 2MM and 2MM + 1MM are also allowed
	--output_dir		path to directory where the output files will be saved


	"""
}

if (params.help) {
	helpMessage()
	exit 0
}

// check if all parameters are provided correctly

if (!file(params.input_bed).exists()) {
	exit 1, "Input bed file not found: ${params.input_bed}"}

if (!file(params.output_dir).exists()) {
	exit 1, "Output directory not found: ${params.output_dir}"}

if (!file(index_bowtie).exists()) {
	exit 1, "Index file not found: ${index_bowtie}"}

if (!file(params.primer_settings).exists()) {
	exit 1, "Primer3 settings file not found: ${params.primer_settings}"}


if (!params.primer3_diff.toString().isNumber()){
	exit 1, "Invalid primer3_diff PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE: ${params.primer3_diff}. Valid options: any integer > 0."}

if (params.primer3_diff.toInteger() < 0){
	exit 1, "Invalid primer3_diff PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE: ${params.primer3_diff}. Valid options: any integer > 0."}

if (!params.primer3_nr.toString().isNumber()){
	exit 1, "Invalid primer3_nr PRIMER_NUM_RETURN: ${params.primer3_nr}. Valid options: any integer > 0. Caution: setting this parameter to a large value will increase running time."}
if (params.primer3_nr.toInteger() < 0){
	exit 1, "Invalid primer3_nr PRIMER_NUM_RETURN: ${params.primer3_nr}. Valid options: any integer > 0. Caution: setting this parameter to a large value will increase running time."}

if (params.snp_filter != "strict" && params.snp_filter != 'loose' && params.snp_filter != 'off'){
	exit 1, "Invalid SNP filter: ${params.snp_filter}. Valid options: 'strict','loose' or 'off'."}

if (params.temp_str_filter != "on" && params.temp_str_filter != 'off'){
	exit 1, "Invalid temp_str_filter filter: ${params.sec_str_filter}. Valid options: 'on','off'."}

if (params.spec_filter != "strict" && params.spec_filter != 'loose'){
	exit 1, "Invalid specificity filter: ${params.spec_filter}. Valid options: 'strict','loose'."}


if (!params.min_tm.toString().isNumber()) {exit 1, " min_tm: ${params.min_tm}. Valid options: any integer > 0."}
if (!params.max_tm.toString().isNumber()) {exit 1, " max_tm: ${params.max_tm}. Valid options: any integer > 0."}
if (!params.opt_tm.toString().isNumber()) {exit 1, " opt_tm: ${params.opt_tm}. Valid options: any integer > 0."}
if (!params.diff_tm.toString().isNumber()) {exit 1, " diff_tm: ${params.diff_tm}. Valid options: any integer > 0."}
if (!params.min_gc.toString().isNumber()) {exit 1, " min_gc: ${params.min_gc}. Valid options: any integer > 0."}
if (!params.max_gc.toString().isNumber()) {exit 1, " max_gc: ${params.max_gc}. Valid options: any integer > 0."}
if (!params.opt_gc.toString().isNumber()) {exit 1, " opt_gc: ${params.opt_gc}. Valid options: any integer > 0."}
if (!params.amp_min.toString().isNumber()) {exit 1, " amp_min: ${params.amp_min}. Valid options: any integer > 0."}
if (!params.amp_max.toString().isNumber()) {exit 1, " amp_max: ${params.amp_max}. Valid options: any integer > 0."}

if (params.min_tm.toInteger() > params.max_tm.toInteger() ) {exit 1, " min_tm and max_tm: max_tm (${params.min_tm}) should be > min_tm (${params.max_tm})"}
if (params.opt_tm.toInteger() > params.max_tm.toInteger() || params.opt_tm.toInteger() < params.min_tm.toInteger() ) {exit 1, " opt_tm: ${params.opt_tm} should > min_tm (${params.min_tm}) and < max_tm (${params.max_tm})"}
if (params.min_gc.toInteger() > params.max_gc.toInteger() ) {exit 1, " min_gc and max_gc: max_gc (${params.min_gc}) should be > min_gc(${params.max_gc})"}
if (params.opt_gc.toInteger() > params.max_gc.toInteger() || params.opt_gc.toInteger() < params.min_gc.toInteger() ) {exit 1, " opt_gc: ${params.opt_gc} should > min_gc (${params.min_gc}) and < max_gc (${params.max_gc})"}


// print run info

log.info """\
==============================================
CIRCprimerXL pipeline
==============================================
OncoRNALab - Marieke Vromman
https://github.com/OncoRNALab/CIRCprimerXL
https://hub.docker.com/repository/docker/oncornalab/CIRCprimerXL
==============================================
your input file: ${params.input_bed}
your output directory: ${params.output_dir}
"""

// define and run each process

process split_circRNAs {

	input:
	path 'input_bed_handle' from input_bed

	output:
	path 'circ*' into ind_circ_file
	path 'start_time.txt' into start_time
	path 'all_circ.txt' into all_cic

	"""
	01_split_input.py -i $input_bed_handle
	python3 -c 'from datetime import datetime; print(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))' > start_time.txt
	"""

}


process get_seq {
	input:
	file ind_circ_file_handle from ind_circ_file.flatten()

	output:
	tuple val("${ind_circ_file_handle.baseName}"), path('input_primer3*') into in_primer3
	tuple val("${ind_circ_file_handle.baseName}"), path('input_NUPACK*') into in_folding
	// tuple val("${ind_circ_file_handle.baseName}"), path('input_SNP*') into in_SNP
	tuple val("${ind_circ_file_handle.baseName}"), path('input_filter*') into in_filter
/* 	tuple val("${ind_circ_file_handle.baseName}"), path('fasta_in.txt') into fasta_SNP
	tuple val("${ind_circ_file_handle.baseName}"), path('fasta_track.txt') into fasta_track
	tuple val("${ind_circ_file_handle.baseName}"), path('annotation*.txt') into annotation_splice */

	"""
	02_get_seq_fastahack.py -i $ind_circ_file_handle -n $params.primer3_diff -p $params.primer3_nr -a $params.min_tm -b $params.max_tm -c $params.opt_tm -d $params.diff_tm -e $params.min_gc -f $params.max_gc -g $params.opt_gc -j $params.amp_min -k $params.amp_max
	"""
}


/* process get_SNPs {

	input:
	tuple val(snp_id), path('in_SNP_handle') from in_SNP
	tuple val(snp_id), path('fasta_SNP_handle') from fasta_SNP
	tuple val(snp_id), path('fasta_track_handle') from fasta_track


	output:
	tuple val(snp_id), path('output_SNP*') into out_SNP_upfront_filter
	tuple val(snp_id), path('output_SNP*') into out_SNP


	"""
	get_SNPs.py -i in_SNP_handle -u $params.snp_url -f fasta_SNP_handle -t fasta_track_handle
	"""
	
} */

process folding_template {

	input:
	tuple val(fold_id_t), path('in_folding_handle') from in_folding

	output:
	tuple val (fold_id_t), path('output_NUPACK_*') into out_folding_template_upfront_filter
	tuple val (fold_id_t), path('output_NUPACK_*') into out_folding_template

	"""
	03_get_sec_str_temp.py -i in_folding_handle -f $params.temp_str_filter
	"""
}



process get_primers {

	publishDir "$params.output_dir/primer3_details", mode: 'copy', pattern: 'output_primer3_*'
	errorStrategy "ignore"

	all_info = out_folding_template_upfront_filter.join(in_primer3).groupTuple()
	
	input:
	tuple val(u_filter_id), path('out_folding_template_upfront_filter_handle'), path('in_primer3_handle') from all_info
/* 	tuple val(u_filter_id), path('out_folding_template_upfront_filter_handle'), path('in_primer3_handle'), path('out_SNP_upfront_filter_handle') from all_info */
	path 'primer_settings_handle' from params.primer_settings
 
	output:
	tuple val(u_filter_id), path('amplicon_folding_input_circ*') into amplicon_folding_in
	tuple val(u_filter_id), path('all_primers_circ*') into all_primers_per_circ
	path('primer_spec_input_circ*') into primer_spec_input_per_circ
	path('output_primer3_*')

	"""
	04_upfront_filter.py -i in_primer3_handle -a out_folding_template_upfront_filter_handle -f $params.temp_str_filter
	/bin/primer3-2.5.0/src/primer3_core --output=output_primer3_${u_filter_id}.txt --p3_settings_file=$primer_settings_handle primer3_file*
	05_split_primers.py -i output_primer3_${u_filter_id}.txt
	"""
}

process folding_amplicon {

	input:
	tuple val(fold_id_a), path(amplicon_folding_in_handle) from amplicon_folding_in.transpose()

	output:
	tuple val(fold_id_a), path('output_NUPACK_*') into out_folding_amplicon

	"""
	06_get_sec_str_amp.py -i $amplicon_folding_in_handle
	"""
}



process specificity_primers {

	// cpus 4
	// cpus 16

	input:
	path 'bowtie_index' from index_bowtie
	file 'primer_spec_input_*' from primer_spec_input_per_circ.collect()

	output:
	file 'out_spec_primer.txt' into out_spec_primer

	script:
	"""
	cat primer_spec_input_* > primer_spec_input_all_circ.txt
	/bin/bowtie-1.3.0-linux-x86_64/bowtie --tryhard -X1000 -v3 --quiet -x bowtie_index/$params.index_bowtie_name --threads ${task.cpus} --12 primer_spec_input_all_circ.txt > out_spec_primer.txt
	"""
}


process filter_primers {

	gather_filter = in_filter.join(out_folding_template).join(out_folding_amplicon).join(all_primers_per_circ).groupTuple()

	input:
	tuple val(filter_id), path('circ_file_handle'), path('out_folding_template_handle'), path('out_folding_amplicon_handle'), path('all_primers_per_circ_handle') from gather_filter
	file 'out_spec_primer_handle' from out_spec_primer

	// val 'snp_filter_handle' from params.snp_filter
	
	output:
	path('selected_primers_*') into results_per_circ
	path('all_primers') into out_dir
	path('log_file_*') into log_file_per_circ

	"""
	mkdir all_primers
	07_filter.py -A circ_file_handle -P all_primers_per_circ_handle -b out_spec_primer_handle -t out_folding_template_handle -a out_folding_amplicon_handle -p $params.spec_filter
	08_gather_output.py -i all_primers/filtered_primers_*
	"""
}



process print_output {

	publishDir params.output_dir, mode: 'copy'

	input:
	file 'results_per_circ*' from results_per_circ.collect()
	file 'log_file_per_circ*' from log_file_per_circ.collect()
	file 'start_time_file' from start_time
	path 'all_primer_files' from out_dir.collect()
	path 'all_circ_file' from all_cic


	output:
	path 'all_primers'
	path 'filtered_primers.txt'
	path 'log_file.txt'
	path 'summary_run.txt'

	"""
	mkdir all_primers
	cp all_primer_files*/* all_primers/
	echo "fusions_seq	seq_id	primer_ID	FWD_primer	REV_primer	FWD_pos	FWD_length	REV_pos	REV_length	FWD_Tm	REV_Tm	FWD_GC	REV_GC	amplicon	PASS" > filtered_primers.txt
	cat results_per_circ* >> filtered_primers.txt
	echo "circ_ID	design	primer_found	total_primer_pairs	passed	failed_spec	failed_sec_str_amp" > log_file.txt
	cat log_file_per_circ* >> log_file.txt
	09_summary_run.py -l log_file.txt -s start_time_file -o . -a all_circ_file
	"""
}
