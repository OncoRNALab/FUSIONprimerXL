# FUSIONprimerXL
Collaborators: Marieke Vromman, Pieter-Jan Volders

Questions concerning the GitHub structure/scripts can be addressed to any of the collaborators.

Primer design pipeline for fusion sequences based on [CIRCprimerXL](https://github.com/OncoRNALab/CIRCprimerXL) (M. Vromman, J. Anckaert, J. Vandesompele, P.J. Volders, CIRCprimerXL: convenient and high-throughput PCR primer design for circular RNA quantification. Frontiers in Bioinformatics (2022), DOI:10.3389/fbinf.2022.834655).

This pipeline runs entirly in the [oncornalab/primerxl_circ](https://hub.docker.com/repository/docker/oncornalab/primerxl_circ) docker image, which is available on DockerHub. It is not necessary to download this image locally, as Nextflow pulls the latest version automatically from DockerHub.

## Installation
### Required reference genome
A reference is required to run the pipeline:
a Bowtie cDNA + ncRNA reference to test the specificity of the primers (https://github.com/BenLangmead/bowtie)

#### Bowtie
In general, a combination of the cDNA and ncRNA index are used to test the specificity of the primers.

```bash
wget ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
gunzip Homo_sapiens.GRCh38.cdna.all.fa.gz
wget ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz
gunzip Homo_sapiens.GRCh38.ncrna.fa.gz
cat Homo_sapiens.GRCh38.cdna.all.fa Homo_sapiens.GRCh38.ncrna.fa > hg38_cdna.fa
rm Homo_sapiens.GRCh38.cdna.all.fa
rm Homo_sapiens.GRCh38.ncrna.fa
bowtie-build hg38_cdna.fa hg38_cdna
```

You do not need to install Bowtie locally to create the required indexes. Instead you can use the Docker container. For this, [Docker](https://docs.docker.com/get-docker/) needs to be installed on your computer.
```bash
docker run -v "$PWD/assets":/assets oncornalab/primerxl_circ:v0.27 /bin/bowtie-1.3.0-linux-x86_64/bowtie-build /assets/index_bowtie/hg38_cdna.fa /assets/index_bowtie/hg38_cdna
```


### Running on your computer
[Nextflow](https://www.nextflow.io/) and [Docker](https://docs.docker.com/get-docker/) should be installed locally. Make sure [Docker Desktop](https://www.docker.com/products/docker-desktop) is running when you want run the pipeline.

### Running on the HPC (UGent)
Nextflow version 20.10.0 is available on all clusters (swalot, skitty, victini, joltik, kirlia, doduo). The pipeline can be run through an interactive session. The pipeline can only run from the $VSC_SCRATCH_VO_USER directory.

```
qsub -I -l nodes=1:ppn=16 -l walltime=04:00:00
cd $VSC_SCRATCH_VO_USER/FUSIONprimerXL/
module load Nextflow/20.10.0
nextflow run FUSIONprimerXL.nf --help
```


## General usage

```
$ nextflow run FUSIONprimerXL.nf --help

Usage:

	The typical command for running the pipeline is as follows:
	nextflow run FUSIONprimerXL.nf -profile singularity

	Mandatory nextflow arguments:
	-profile 		set to 'local' when running locally, set to 'singularity' when running on the HPC

	Mandatory pipeline arguments:
	--input_seq			path to input file with fusion sequences (the fusion should be indicated with a *)
	--index_bowtie		path to bowtie genome index directory
	--index_bowtie_name	the basename of the Bowtie index to be searched (the name of any of the index files up to but not including the final .1.ebwt / .rev.1.ebwt / ...)


	Optional pipeline arguments:
	--primer_settings	path to file with primer3plus settings (see primer3 manual)
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
	--temp_str_filter	when set to 'on', secunday structures in the tample are avoided during primer design; when set to 'off', this is not done (default: 'on')
	--spec_filter		when set to 'strict', only 2MM + 3MM are allowed; when set to 'loose', 2MM + 2MM and 2MM + 1MM are also allowed
	--output_dir		path to directory where the output files will be saved

	
```

You can easily create your own profiles by modifying the nextflow.config file.

Nextflow keeps track of all the processes executed in your pipeline. If you want to restart the pipeline after a bug, add '-resume'. The execution of the processes that are not changed will be skipped and the cached result used instead.

Note: The pipeline results are cached by default in the directory $PWD/work. This folder can take of lot of disk space. If your are sure you won’t resume your pipeline execution, clean this folder periodically.

## Output
In the output folder, you will find
<ul>
  <li>filtered_primers.txt, a file containing one selected primer pair per fusion (see below for column details)</li>
  <li>log_file.txt, </li>
  <li>summary_run.txt </li>
  <li>all_primers directory</li>
  <li>primer3_details directory</li>
</ul>

filtered_primers.txt output file column names:

| column name      | description                                                                                                            |
|:-----------------|:-----------------------------------------------------------------------------------------------------------------------|
| fusion_ID          | fusion id assigned to each fusion sequence (unique within one run)                                                               |
| seq_id              | fusion sequence                                                                                                        |
| primer_ID        | primer ID generated by primer3                                                                                         |
| FWD_primer       | forward primer                                                                                                         |
| REV_primer       | reverse primer                                                                                                         |
| FWD_pos          | relative position of forward primer                                                                                    |
| FWD_length       | length of forward primer                                                                                               |
| REV_pos          | relative position of reverse primer                                                                                    |
| REV_length       | length of reverse primer                                                                                               |
| FWD_Tm           | melt temperature of forward primer                                                                                     |
| REV_Tm           | melt temperature of reverse primer                                                                                     |
| FWD_GC           | GC content of forward primer                                                                                           |
| REV_GC           | GC content of reverse primer                                                                                           |
| amplicon         | amplicon sequence amplified by the primer pair                                                                         |
| PASS             | result of filtering (PASS if the primer pair passed all filters, FAIL if   the primer pair failed one or more filters) |


## Other species
As default, FUSIONprimerXL designs primers voor humans. To design primers for other species, the following files have to be provided and parsed through the corresponding parameters:
<ul>
  <li>a file containing the chromosome sizes (parameter chrom_file) (for example from: https://www.ncbi.nlm.nih.gov/grc/human/data)</li>
  <li>the correct bowtie index mentioned above
</ul>


## Nextflow tower

[Nextflow tower](https://tower.nf/) can be used to monitor the pipeline while it's running.
```
nextflow run FUSIONprimerXL.nf -with-tower
```

When Nextflow tower is used in combination with the HPC, the nextflow version and tower access token should be indicated.
```
export NXF_VER="20.10.0"
export TOWER_ACCESS_TOKEN=your_token_here
```

