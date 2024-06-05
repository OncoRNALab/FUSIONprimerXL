#!/usr/bin/env nextflow

// Enable DSL2
nextflow.enable.dsl=2

/*
====================================================================================================
Pipeline: FUSIONprimerXL
Description: A Nextflow pipeline for the design of fusion transcript-specific primers. 
License: MIT
Copyright (c) 2021 Ghent University
====================================================================================================
*/
/*
====================================================================================================
DEFAULT PARAMETERS (can be overwrittien in config file)
====================================================================================================
*/
// Input file
params.input_bed = "$projectDir/input/input_fusionRNAs.bed"
// Prediction program
params.prediction_program = "ViennaRNA"
// Bowtie index
params.index_bowtie = "$projectDir/assets/GRCh38/index_bowtie"
params.index_bowtie_name = "hg38_cdna"
// Fastahack index
params.index_fasta = "$projectDir/assets/GRCh38/index_fastahack"
params.index_fasta_name = "Homo_sapiens.GRCh38.dna.primary_assembly.fa"
// Primer3 settings
params.primer_settings = "$projectDir/assets/primer3_settings.txt"
// chromosome sizes
params.chrom_file = "$projectDir/assets/GRCh38/chrom_sizes_GRCh38.txt"
// known exons
params.known_exons = "$projectDir/assets/GRCh38/Known_exons_GRCh38.111.bed"
// ENST_list
params.list_ENST = "$projectDir/assets/GRCh38/ENST_list_GRCh38.txt"
// output directory
params.output_dir = "$projectDir/output/"
// splice
params.splice = 'yes'
// Primer3 parameters:
params.primer3_diff = 1
params.primer3_nr = 20
// Primer Tm parameters:
params.min_tm = 58
params.opt_tm = 59
params.max_tm = 60
params.diff_tm = 2
// Primer GC parameters:
params.min_gc = 30
params.opt_gc = 50
params.max_gc = 80
// Amplicon length parameters:
params.amp_min = 50
params.amp_max = 0 // this param is set to 0, so that it can be adjusted depending on temp_l if the user does not supply amp_max
params.temp_l = 150

// filters
params.snp_filter = 'strict'
params.spec_filter = 'strict'
params.upfront_filter = "yes"
// SNP database
params.snp_url = 'http://hgdownload.soe.ucsc.edu/gbdb/hg38/snp/dbSnp155Common.bb'

// display help message
params.help = false

// required parameters
input_bed = file(params.input_bed)
chrom_file = file(params.chrom_file)
index_bowtie = file(params.index_bowtie)
index_fasta = file(params.index_fasta)
known_exons = file(params.known_exons)
list_ENST = file(params.list_ENST)


/*
====================================================================================================
HELP MESSAGE
====================================================================================================
*/

def helpMessage() {
	log.info"""
	
	Usage:
	
	The typical command for running the pipeline is as follows (standard = default parameters):
	nextflow run FUSIONprimerXL.nf -profile standard 
	
	Mandatory nextflow arguments:
	-profile 		set to 'local' when running locally, set to 'singularity' when running on the HPC

	
	Mandatory pipeline arguments:
	--input_bed			path to input file with fusionRNAs in bed format (0-based annotation)
	--index_bowtie		path to bowtie genome index directory
	--index_bowtie_name	the basename of the Bowtie index to be searched (the name of any of the index files up to but not including the final .1.ebwt / .rev.1.ebwt / ...)
	--index_fasta		path to directory that contains the fastahack genome and index file
	--index_fasta_name	the name of the fastahack genome file

	Optional pipeline arguments:
	--prediction_program	program used for secondary structure prediction [ViennaRNA or Nupack](default: ViennaRNA)
	--splice			when set to 'yes' the input sequence will be spliced, when set to 'no' the input sequence will be unspliced
	--primer_settings	path to file with primer3 settings (see primer3 manual)
	--chrom_file		file containing all chromosome sizes (to validate bed file)
	--known_exons		bed file containing exon annotation
	--list_ENST			file containing ENST numbers of canonical transcripts or transcripts of interest (this file can also be left empty)
	--primer3_diff		the minimum number of base pairs between the 3' ends of any two left primers (see also primer3 PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE)
	--primer3_nr		the number of primers designed by primer3; caution: setting this parameter to a large value will increase running time
	--min_tm			minimum melt temperature of the primers (default: 58)
	--max_tm			maximum melt temperature of the primers(default: 60)
	--opt_tm			optimal melt temperature of the primers(default: 59)
	--diff_tm			maximum difference in melt temperature between the primers(default: 2)
	--min_gc			minimum GC contect of the primers (default: 30)
	--max_gc			maximum GC contect of the primers(default: 80)
	--opt_gc			optimal GC contect of the primers(default: 50)
	--amp_min			minimum amplicon length (default: 60)
	--amp_max			maximum amplicon length (default: 0)
	--temp_l			the number of nucleotides on each side of the breakpoint that will be used for the template (example 150 => template of 300 nts in total)
	--spec_filter		when set to 'strict', a maximum of 2 MM is allowed; when set to 'loose', a maximum of 5 combined MM is allowed
	--snp_filter		when set to 'strict', no common SNPs are allowed in primer sequence; when set to 'loose', common SNPs are allowed in 5' first half of primer; when set to 'off', no filter is applied
	--snp_url			when using a differente species than human, the correct SNP database url should be provided; alternatively, this paramater can be set to 'off' if no SNP database is available
	--upfront_filter	when set to 'yes', SNPs and secundary structures are avoided before primer design; when set to 'str', secundary structures are avoided before primer design; when set to 'snp', snp are avoided before primer design; when set to 'no', no filtering before primer design is performed
	--output_dir		path to directory where the output files will be saved


	"""
}
// help message
if (params.help) {
	helpMessage()
	exit 0
}

/*
====================================================================================================
CHECK IF ALL PARAMETERS ARE PROVIDED CORRECTLY
====================================================================================================
*/
// --input_seq
if (!file(params.input_bed).exists()) {exit 1, "Input bed file not found: ${params.input_bed}"}
// --out_dir
if (!file(params.output_dir).exists()) {exit 1, "Output directory not found: ${params.output_dir}"}
// --index_bowtie
if (!file(index_bowtie).exists()) {exit 1, "Index file not found: ${index_bowtie}"}
// --primer_settings
if (!file(params.primer_settings).exists()) {exit 1, "Primer3 settings file not found: ${params.primer_settings}"}
// --chrom_file
if (!file(chrom_file).exists()) {exit 1, "Chromosome file not found: ${chrom_file}"}
// -- splice
if (params.splice != "yes" && params.splice != 'no'){
	exit 1, "Invalid splicing option: ${params.splice}. Valid options: 'yes','no'."}
if (params.splice == "yes"){if (!file(params.known_exons).exists()) {exit 1, "Known exons file not found: ${params.known_exons}"}}
if (params.splice == "yes" && !(params.list_ENST == "none")) {if (!file(params.list_ENST).exists()) {exit 1, "ENST list file not found: ${params.list_ENST}"}}
// --primer3_diff
if (!params.primer3_diff.toString().isNumber()){
	exit 1, "Invalid primer3_diff PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE: ${params.primer3_diff}. Valid options: any integer > 0."}
if (params.primer3_diff.toInteger() < 0){
	exit 1, "Invalid primer3_diff PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE: ${params.primer3_diff}. Valid options: any integer > 0."}
// --primer3_nr
if (!params.primer3_nr.toString().isNumber()){
	exit 1, "Invalid primer3_nr PRIMER_NUM_RETURN: ${params.primer3_nr}. Valid options: any integer > 0. Caution: setting this parameter to a large value will increase running time."}
if (params.primer3_nr.toInteger() < 0){
	exit 1, "Invalid primer3_nr PRIMER_NUM_RETURN: ${params.primer3_nr}. Valid options: any integer > 0. Caution: setting this parameter to a large value will increase running time."}
// --temp_l
if (!params.temp_l.toString().isNumber()){
	exit 1, "Invalid template length: ${params.temp_l}. Valid options: any integer between 50 and 500."}
if (params.temp_l.toInteger() < 50 || params.temp_l.toInteger() > 500){
	exit 1, "Invalid template length: ${params.temp_l}. Valid options: any integer between 50 and 500."}
// --snp_filter
if (params.snp_filter != "strict" && params.snp_filter != 'loose' && params.snp_filter != 'off'){
	exit 1, "Invalid SNP filter: ${params.snp_filter}. Valid options: 'strict','loose'."}
// --spec_filter
if (params.spec_filter != "strict" && params.spec_filter != 'loose'){
	exit 1, "Invalid specificity filter: ${params.spec_filter}. Valid options: 'strict','loose'."}
// --upfront_filter
if (params.upfront_filter != "yes" && params.upfront_filter != 'str' && params.upfront_filter != 'snp' && params.upfront_filter != 'no'){
	exit 1, "Invalid SNP filter: ${params.upfront_filter}. Valid options: 'yes','str','snp','no'."}
//	--min_tm
if (!params.min_tm.toString().isNumber()) {exit 1, " min_tm: ${params.min_tm}. Valid options: any integer > 0."}
//	--max_tm
if (!params.max_tm.toString().isNumber()) {exit 1, " max_tm: ${params.max_tm}. Valid options: any integer > 0."}
//	--opt_tm
if (!params.opt_tm.toString().isNumber()) {exit 1, " opt_tm: ${params.opt_tm}. Valid options: any integer > 0."}
//	--diff_tm
if (!params.diff_tm.toString().isNumber()) {exit 1, " diff_tm: ${params.diff_tm}. Valid options: any integer > 0."}
//	--min_gc
if (!params.min_gc.toString().isNumber()) {exit 1, " min_gc: ${params.min_gc}. Valid options: any integer > 0."}
//	--max_gc
if (!params.max_gc.toString().isNumber()) {exit 1, " max_gc: ${params.max_gc}. Valid options: any integer > 0."}
//	--opt_gc
if (!params.opt_gc.toString().isNumber()) {exit 1, " opt_gc: ${params.opt_gc}. Valid options: any integer > 0."}
//	--amp_min
if (!params.amp_min.toString().isNumber()) {exit 1, " amp_min: ${params.amp_min}. Valid options: any integer > 0."}
//	--amp_max
if (!params.amp_max.toString().isNumber()) {exit 1, " amp_max: ${params.amp_max}. Valid options: any integer >= 0, input 0 to adjust depending on temp_l."}

//	checking logic
if (params.min_tm.toInteger() > params.max_tm.toInteger() ) {exit 1, " min_tm and max_tm: max_tm (${params.min_tm}) should be > min_tm (${params.max_tm})"}
if (params.opt_tm.toInteger() > params.max_tm.toInteger() || params.opt_tm.toInteger() < params.min_tm.toInteger() ) {exit 1, " opt_tm: ${params.opt_tm} should > min_tm (${params.min_tm}) and < max_tm (${params.max_tm})"}
if (params.min_gc.toInteger() > params.max_gc.toInteger() ) {exit 1, " min_gc and max_gc: max_gc (${params.min_gc}) should be > min_gc(${params.max_gc})"}
if (params.opt_gc.toInteger() > params.max_gc.toInteger() || params.opt_gc.toInteger() < params.min_gc.toInteger() ) {exit 1, " opt_gc: ${params.opt_gc} should > min_gc (${params.min_gc}) and < max_gc (${params.max_gc})"}

/*
====================================================================================================
RUN INFO
====================================================================================================
*/
log.info """\
==============================================
FUSIONprimerXL pipeline
==============================================
OncoRNALab - Arne Blom / Marieke Vromman 
Github -
Docker - 
==============================================
your input file: ${params.input_bed}
your output directory: ${params.output_dir}
"""
/*
====================================================================================================
PROCESS 1 - splitting input bed file with fusionRNAs
====================================================================================================
*/
// channels
input_bed_handle = channel.fromPath(params.input_bed)
chrom_file_handle = channel.fromPath(params.chrom_file)

// process
process split_fusionRNAs {
	tag "split_fusionRNAs"
	input:
	path('input_bed_handle')
	path('chrom_file_handle')

	output:
	path 'fusion*'
	path 'start_time.txt'
	path 'all_fusions.txt'

	"""
	00_validate_bed.py -i $input_bed_handle -c $chrom_file_handle
	01_split_fusionRNAs.py -i $input_bed_handle
	python3 -c 'from datetime import datetime; print(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))' > start_time.txt
	"""
}

/*
====================================================================================================
PROCESS 2 - get_sequence of fusionRNAs
====================================================================================================
*/

// process
process get_sequence {
	tag "get_sequence"

	input:
	path(ind_fusion_file_handle)
	path('known_exons')
	path('list_ENST')
	path('index_fasta')

	output:
	tuple val("${ind_fusion_file_handle.getBaseName()}"), path('annotation*.txt')
	tuple val("${ind_fusion_file_handle.getBaseName()}"), path('fasta_track.txt')
	tuple val("${ind_fusion_file_handle.getBaseName()}"), path('fasta_in.txt')
	tuple val("${ind_fusion_file_handle.getBaseName()}"), path('input_filter*')
	tuple val("${ind_fusion_file_handle.getBaseName()}"), path('input_get_sec_structure_*')
	tuple val("${ind_fusion_file_handle.getBaseName()}"), path('input_primer3_*')
	tuple val("${ind_fusion_file_handle.getBaseName()}"), path('input_SNP_*')
	tuple val("${ind_fusion_file_handle.getBaseName()}"), path('fasta_out.txt')

	
	"""
	02_Sequence_retrieval_chrom_pos.py -n $params.temp_l -i $ind_fusion_file_handle -s $params.splice -e $known_exons -t $list_ENST
	cat fasta_in.txt | fastahack -c $index_fasta/$params.index_fasta_name > fasta_out.txt
	03_generate_in_files.py -i $ind_fusion_file_handle -p $params.primer3_nr -n $params.primer3_diff -a $params.min_tm -b $params.max_tm -c $params.opt_tm -d $params.diff_tm -e $params.min_gc -f $params.max_gc -g $params.opt_gc -j $params.amp_min -k $params.amp_max
	"""
}

/*
====================================================================================================
PROCESS 3 - get single nucleotide polymorphisms
====================================================================================================
*/

// process
process get_SNPs {
	tag "get_SNPs"

	input:
	tuple val(snp_id), path('in_SNP_handle')
	tuple val(snp_id), path('fasta_SNP_handle')
	tuple val(snp_id), path('fasta_track_handle')

	output:
	tuple val(snp_id), path('output_SNP_*')

	"""
	04_get_SNPs.py -i $in_SNP_handle -u $params.snp_url -f $fasta_SNP_handle -t $fasta_track_handle
	"""
}

/*
====================================================================================================
PROCESS 4 - folding_template
====================================================================================================
*/
mode = params.prediction_program

process folding_template {
	tag "folding_template"

	input:
	tuple val(fold_id_t), path('inp_fold_handle')

	output:
	tuple val(fold_id_t), path('output_RNAfold_temp_*')

	script:
	if ( mode == "ViennaRNA") 
		"""
		05_get_sec_str_temp.py -i $inp_fold_handle
		"""
	else if ( mode == "Nupack")
		"""
		Nupack_temp.py -i $inp_fold_handle
		"""
}

/*
====================================================================================================
PROCESS 5 - upfront filter
====================================================================================================
*/

process upfront_filter{
	tag "upfront_filter"

	input:
	tuple val(upfront_filter_id), path(input_primer3_handle), path(out_RNAfold_temp_handle), path(input_SNP_handle)

	output:
	tuple val(upfront_filter_id), path('primer3_file_*')

	"""
	06_upfront_filter.py -i $input_primer3_handle -a $out_RNAfold_temp_handle -s $input_SNP_handle -f $params.upfront_filter -l $params.temp_l
	"""
}

/*
====================================================================================================
PROCESS 6 - primer3
====================================================================================================
*/

process get_primers {
	tag "get_primers"
	publishDir "$params.output_dir/primer3_details", mode: 'copy', pattern: 'output_primer3_*'

	input:
	tuple val(primer3_id), path(primer3_file_handle)
	path(primer_settings_handle)

	output:
	tuple val(primer3_id), path('output_primer3_*')

	
	"""
	06_primer3.py -i $primer3_file_handle -x $primer_settings_handle
	"""
}

/*
====================================================================================================
PROCESS 7 - split primers
====================================================================================================
*/
process split_primers {
	tag "split_primers"
	
	input:
	tuple val(primer3_id), path(output_primer3_handle)

	output:
	tuple val(primer3_id), path('amplicon_folding_input_fusion*')
	tuple val(primer3_id), path('all_primers_fusion*')
	tuple val(primer3_id), path('primer_spec_input_fusion*')

	"""
	07_split_primers.py -i $output_primer3_handle
	"""
}

/*
====================================================================================================
PROCESS 8 - folding_amplicon
====================================================================================================
*/

process folding_amplicon {
	tag "folding_amplicon"

	input:
	tuple val(fold_amp_id), path(amplicon_folding_input_handle)

	output:
	tuple val(fold_amp_id), path('output_RNAfold_amp_*')

	script:
	if ( mode == "ViennaRNA") 
		"""
		08_get_sec_str_amp.py -i $amplicon_folding_input_handle
		"""

	else if ( mode == "Nupack")
		"""
		Nupack_amp.py -i $amplicon_folding_input_handle
		"""
}

/*
====================================================================================================
PROCESS 9 - primer specificity
====================================================================================================
*/

process specificity_primers {
	tag "specificity_primers"

	input:
	tuple val(spec_id), path(primer_spec_input_handle)
	path(bowtie_index)

	output:
	tuple val(spec_id), path('primer_spec_input_all_fusion.txt')
	tuple val(spec_id), path('out_spec_primer.txt')

	"""
	cat $primer_spec_input_handle > primer_spec_input_all_fusion.txt
	bowtie --tryhard -X1000 -v3 --quiet -x $bowtie_index/$params.index_bowtie_name --threads ${task.cpus} --12 primer_spec_input_all_fusion.txt > out_spec_primer.txt
	"""
}

/*
====================================================================================================
PROCESS 10 - filter_primers
====================================================================================================
*/

process filter_primers {
	tag "filter_primers"

	input:

		tuple val(filter_id), path(input_filter_handle), path(all_primers_handle), path(output_SNP_handle), path(output_RNAfold_temp_handle), path(output_RNAfold_amp_handle), path(out_spec_primer_handle), path(annotation_splice_handle)

	output:

	path('selected_primers_*')
	path('all_primers')
	path('log_file_*')

	"""
	mkdir all_primers
	09_filter.py -A $input_filter_handle -P $all_primers_handle -l 150 -s $output_SNP_handle -t $output_RNAfold_temp_handle -a $output_RNAfold_amp_handle -b $out_spec_primer_handle -p $params.spec_filter -f $params.snp_filter
	10_gather_output.py -i all_primers/all_* -a $annotation_splice_handle
	"""
}

/*
====================================================================================================
PROCESS 11 - print_output
====================================================================================================
*/

process print_output {

	tag "print_output"
	publishDir params.output_dir, mode: 'copy'

	input:
	path('results_per_fusion*')
	path('log_file_per_fusion*')
	path('start_time_file')
	path('all_primer_files')
	path('all_fusions_file')


	output:
	path('all_primers')
	path('suggested_primer_pairs.txt')
	path('log_file.txt')
	path('summary_run.txt')

	"""
	mkdir all_primers
	cp all_primer_files*/* all_primers/
	echo "fusion_ID\tprimer_ID\tFWD\tREV\tFWD_posistion\tFWD_length\tREV_position\tREV_length\tprimer_left_TM\tprimer_right_TM\tprimer_left_GC_perc\tprimer_right_GC_perc\tamplicon\tfilter\tleft_annotation\tright_annotation \tsplicing" > suggested_primer_pairs.txt
	cat results_per_fusion* >> suggested_primer_pairs.txt
	echo "fusion_ID\tchrom1\tend\tchrom2\tstart\tdesign\tprimer_found\ttotal_primers\tpassed\tfailed_spec\tfailed_SNP\tfailed_str_temp\tfailed_str_amp" > log_file.txt
	cat log_file_per_fusion* >> log_file.txt
	11_summary_run.py -l log_file.txt -s start_time_file -o . -u $params.upfront_filter -a all_fusions_file
	"""
}

/*
====================================================================================================
THE WORKFLOW
====================================================================================================
*/
workflow {
	split_fusionRNAs(
		input_bed_handle,
		chrom_file_handle
	)

	get_sequence(
		//ind_fusion_file_handle
		split_fusionRNAs.out[0].flatten(),
		known_exons,
		list_ENST,
		index_fasta
	)

	get_SNPs (
		// in_SNP_handle
		get_sequence.out[6],
		// fasta_SNP_handle
		get_sequence.out[2],
		// fasta_track_handle
		get_sequence.out[1]
	)

	folding_template (
		// inp_fold_handle
		get_sequence.out[4]
	)

	upfront_filter (
		get_sequence.out[5].join(folding_template.out[0]).join(get_SNPs.out[0]).groupTuple()
	)

	get_primers (
		upfront_filter.out[0],
		params.primer_settings
	)

	split_primers (
		get_primers.out[0]
	)

	folding_amplicon (
		split_primers.out[0]
	)

	specificity_primers (
		split_primers.out[2],
		index_bowtie
	)

	filter_primers (
		get_sequence.out[3].join(split_primers.out[1]).join(get_SNPs.out[0]).join(folding_template.out[0]).join(folding_amplicon.out[0]).join(specificity_primers.out[1]).join(get_sequence.out[0]).groupTuple()
	)

	print_output(
		filter_primers.out[0].collect(),
		filter_primers.out[2].collect(),
		split_fusionRNAs.out[1],
		filter_primers.out[1].collect(),
		split_fusionRNAs.out[2]
	)
}

workflow.onComplete {
	println "\n\t\t\t  Pipeline execution summary\n"+
		"=================================================================================\n\n"+
		"\tPipeline completed at:\t$workflow.complete\n" +
		"\tExecution status:\t${ workflow.success ? 'OK' : 'failed' }\n"+
		"\tNextflow version:\t$nextflow.version\n"+
		"\tExecuted command:\t$workflow.commandLine\n"+
		"\tRun name:\t\t$workflow.runName\n"+
		"\tConfig file:\t\t$workflow.configFiles\n"+
		"\tProfile:\t\t$workflow.profile\n"+
		"\tContainer engine:\t$workflow.containerEngine\n"+
		"\tContainer:\t\t$workflow.container\n"+
		"\tStart time:\t\t$workflow.start\n"+
		"\tCompletion:\t\t$workflow.complete\n"+
		"\tDuration:\t\t$workflow.duration\n"+
		"\tProject directory:\t$workflow.projectDir\n"+
		"\tExit status:\t\t$workflow.exitStatus\n"+
		"================================================================================="
}