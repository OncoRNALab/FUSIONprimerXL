# FUSIONprimerXL
<hr>

### Info

- Collaborators: Marieke Vromman, Arne Blom, Annelien Morlion. Questions concerning the GitHub structure/scripts can be addressed to any of the collaborators.
- Primer design pipeline for fusionRNAs based on primerXL (Lefever, S., Pattyn, F., De Wilde, B. et al. High-throughput PCR assay design for targeted resequencing using primerXL. BMC Bioinformatics 18, 400 (2017), https://doi.org/10.1186/s12859-017-1809-3) and CIRCprimerXL (M. Vromman, J. Anckaert, J. Vandesompele, P.J. Volders, CIRCprimerXL: convenient and high-throughput PCR primer design for circular RNA quantification. Frontiers in Bioinformatics (2022), https://doi.org/10.3389/fbinf.2022.834655).
- This pipeline runs entirly in the [oncornalab/fusionprimerxl](https://hub.docker.com/r/oncornalab/fusionprimerxl/tags) docker image, which is available on DockerHub. It is not necessary to download this image locally, as Nextflow pulls the latest version automatically from DockerHub (depending on the profile).
 

### Table of contents

- [1. Installation](#1-Installation)
	- [Fastahack index](#Building-fastahack-index)
	- [Bowtie index](#Building-the-Bowtie-index)
- [Running](#2-Running-on-your-computer)
	- [Example](#Example)
- [HPC](#3-Step-Running-on-the-HPC-UGent)
- [Other species](#4-Other-species)
- [Nextflow tower](#5-Nextflow-tower)
- [CITE](#Cite)

## 1. Installation
<hr>

### Cloning the repository:
```
git clone https://github.com/OncoRNALab/FUSIONprimerXL.git
```

### Requirements:
The pipeline can be run entirely on the  [oncornalab/fusionprimerxl](https://hub.docker.com/r/oncornalab/fusionprimerxl) docker image, which will automatically be pulled by Nextflow if specified in the profile. Two refences are required to run the pipeline, these are not included in the git repository because of their size. Instructions on how to build them can be found in [getting started](#Getting-started). The pipeline can also be tested with the included [Example](#example) (example directory) without having to build any references, afterwards references can be build to make use of the pipeline.

#### A) Local software:

Programs that need to be installed locally:
- [Nextflow](https://www.nextflow.io/docs/latest/install.html)
- [Docker](https://docs.docker.com/engine/install/)

#### B) References:

The following references (2 indexes and 3 files) are required to run the pipeline:
- a fasthack reference genome to extract the fusionRNA sequence surrounding the breakpoint (https://github.com/ekg/fastahack) (**needs to be build see below**)
- a Bowtie cDNA + ncRNA reference to test the specificity of the primers (https://github.com/BenLangmead/bowtie) (**needs to be build see below**)
- a bed file containing all exons. (included for for GRCh38)
- a file containing the chromosome sizes (**included** in assets)
- a file containing the canonical transcripts (**included** in assets/GRCh38)

Both indexes below are examples and can be replaced by your choosen index (--index_fasta, --index_bowtie) and species (see below).

### Getting started:

You do not need to install fastahack and Bowtie locally to create the required indexes. Instead you can use the Docker container. For this, [Docker](https://docs.docker.com/get-docker/) needs to be installed on your computer.

### <u>Building fastahack index</u>

Make sure to to download the fastafile in a folder called index_fastahack for the corresponding genome (example : /assets/GRCh38/index_fastahack/). The steps for GRCh38 are described below, for other species see [other species](#4-Other-species), replace the corresponding links and file locations.
- Step 1:  download the complete primary assembly in the corresponding assets folder.

```
wget https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
```
- Step 2: Unzip the file

```
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
```
- Step 3: Create the index with fastahack (make sure to run this command from the base folder FUSIONprimerXL)

***With the docker image***:
```
docker run -v "$PWD":/assets oncornalab/fusionprimerxl:v0.2 fastahack -i assets/Homo_sapiens.GRCh38.dna.primary_assembly.fa
```

***Or locally***:
```
fastahack -i assets/GRCh38/index_fastahack/Homo_sapiens.GRCh38.dna.primary_assembly.fa
```

- Step4: Move the files to the correct assets folder.
```
mv ./Homo_sapiens.GRCh38.dna.primary_assembly.fa ./assets/GRCh38/index_fastahack/Homo_sapiens.GRCh38.dna.primary_assembly.fa ;
mv ./Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai ./assets/GRCh38/index_fastahack/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai
```

### <u>Building the Bowtie index</u>

In general, a combination of the cDNA and ncRNA index are used to test the specificity of the primers.  The steps for GRCh38 are described below, for other species see [other species](#4-Other-species), replace the corresponding links and file locations (example: /assets/GRCh38/index_bowtie/).

- Step 1: download the cDNA and ncRNA fasta files in the corresponding assets folder.

```
wget https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz ; wget https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz
```

- Step 2: unzip the downloaded fasta file.

```
gunzip Homo_sapiens.GRCh38.cdna.all.fa.gz ; gunzip Homo_sapiens.GRCh38.ncrna.fa.gz
```

- Step 3: write the fasta files to a single file.

```
cat Homo_sapiens.GRCh38.cdna.all.fa Homo_sapiens.GRCh38.ncrna.fa > hg38_cdna.fa ; rm Homo_sapiens.GRCh38.cdna.all.fa ; rm Homo_sapiens.GRCh38.ncrna.fa
```

- Step 4: build the index

***With docker (make sure to run the command from the base folder FUSIONprimerXL)***
```
docker run -v "$PWD":/assets oncornalab/fusionprimerxl:v0.2 bowtie-build /assets/hg38_cdna.fa /assets/hg38_cdna
```

***Or locally***
```
bowtie-build ./assets/GRCh38/index_bowtie/hg38_cdna.fa ./assets/GRCh38/index_bowtie/hg38_cdna
```

- Step 5: move the files to the correct location:
```
mv ./hg38_cdna* ./assets/GRCh38/index_bowtie/
```

### Bed file with canonical exons (included for GRCh38)
This file is included for GRCh38; however, the process is the same for  [other species](#4-Other-species), the steps are described below.

Step 1: download the gtf file from ensemble (https://www.ensembl.org/info/data/ftp/index.html) in the corresponding assets folder (example: assets/GRch38/)
```
wget https://ftp.ensembl.org/pub/release-111/gtf/homo_sapiens/Homo_sapiens.GRCh38.111.gtf.gz; gunzip Homo_sapiens.GRCh38.111.gtf.gz
```

Step 2: transform the gtf file to a bed file (execute from base folder FUSIONprimerXL)
***NOTE: if you are doing this for C. elegans use ./bin/00_A_gtf_to_bed_C_elegans.py instead***
```
python3 ./bin/00_A_gtf_to_bed.py -i ./Homo_sapiens.GRCh38.111.gtf -o ./assets/GRCh38/Known_exons_GRCh38.111.bed ; rm ./Homo_sapiens.GRCh38.111.gtf
```


Step 3: validate the bed file. This will probably show the warning faulty chromosome name; check the file at the specified lines if the chromosome name is indeed incorrect, correct it. If chromosome names like chrMT or chrKI270728.1 are present these lines can be removed or ignored.
```
python3 ./bin/00_B_validate_bed.py -i ./assets/GRCh38/Known_exons_GRCh38.111.bed -c ./assets/GRCh38/chrom_sizes_GRCh38.txt 
```

### A file containing the canonical transcripts (**included** in assets/GRCh38)
This file is already included and does not need to be created, the steps below are described for your information (FYI).  (https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/):
```
# step 1 get the MANE file (in the correct assets folder)
wget https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/MANE.GRCh38.v1.3.ensembl_genomic.gtf.gz ; gunzip MANE.GRCh38.v1.3.ensembl_genomic.gtf.gz
# step 2 generate ENST_list
python3 ./bin/00_C_generate_ENST_list.py -i ./assets/GRCh38/MANE.GRCh38.v1.3.ensembl_genomic.gtf -o ./assets/GRCh38/ENST_list_GRCh38.txt 
```

### The assets folder should now look like this

- assets/GRCh38/
	- index_fastahack/
	- index_bowtie/
	- ENST__list_GRCh38.txt
	- Known_exons_GRCh38.111.bed
	- chrom_sizes_GRCh38.txt

## 2. Running on your computer
<hr>

[Nextflow](https://www.nextflow.io/docs/latest/install.html) and [Docker](https://docs.docker.com/get-docker/) should be installed locally. Make sure the docker engine or [Docker Desktop](https://www.docker.com/products/docker-desktop) is running when you want run the pipeline.

### Input

This tool requires a bed file as input, the input can have a User-defined filename. Good practice would be to place your input files in the input directory. The location of the input file can be defined in the nextflow.config file as params.input_bed = "Path" or by using the --input_bed \<path> on the command line.
The input file with end and start representing the breakpoint:

Note: make sure to use tabs instead of four spaces

Important: make sure to use **0-based counting**

> chromsome1    end    chromsome2    start
```
chr6	160650497	chr6	160666228
chr21	37840709 	chr21	38121128
chr2	44270113 	chr2	237973225
chr5	70901933	chr5	70938841
```

### Output

In the output folder, you will find:
<ul>
  <li>suggested_primer_pairs.txt, a file containing one selected primer pair per fusionRNA (see below for column details)</li>
  <li>log_file.txt, a file containing some statistics per fusionRNA</li>
  <li>summary_run.txt, a file containing general statics for the whole run. </li>
  <li>all_primers directory, contains files with all primers that were created per fusionRNA</li>
  <li>primer3_details directory, contains files with the primer3 details for each fusionRNA</li>
</ul>

filtered_primers.txt output file column names:

| column name      | description                                                                                                            |
| :--------------- | :--------------------------------------------------------------------------------------------------------------------- |
| fusion_ID        | fusion id assigned to each fusionRNA (unique within one run)                                                           |
| primer_ID        | primer ID generated by primer3 (unique per fusion ID)                                                                  |
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
| left_annotation  | exons used for the template sequence on the right side of the BP                                                       |
| right_annotation | exons used for the template sequence on the left side of the BP                                                        |
| splicing         | spliced or unspliced template sequence was used                                                                        |


### Example

This repository contains an example run with 3 fusionRNAs. For this, a small subset of the indexes are also present in the example folder. The example can be run by:
```
nextflow run FUSIONprimerXL.nf -profile example
```

### General usage
From this point onward you are using the full tool, which means you will have to build the corresponding indexes as described in [getting started](#Getting-started). Only the example can be ran without any additional references.

To display information about all available parameters
`nextflow run FUSIONprimerXL.nf --help`

A command can be run like the example below or by defining the parameters in nextflow.config. ***Note: index_bowtie, index_bowtie_name, index_fasta and index_fasta_name have default values for GRCh38 corresponding to the steps above. If you made any changes or are working with another organism make sure the names and locations are correct for these arguments!***
```
nextflow run FUSIONprimerXL.nf -profile standard --input_bed $projectDir/input/input_fusionRNAs.bed --index_bowtie $projectDir/assets/index_bowtie/ --index_bowtie_name hg38_cdna --index_fasta $projectDir/assets/GRCh38/index_fastahack --index_fasta_name Homo_sapiens.GRCh38.dna.primary_assembly.fa
```

All parameters:

```
Usage:

        The typical command for running the pipeline is as follows (standard = default parameters):
        nextflow run FUSIONprimerXL.nf -profile standard 

        Mandatory nextflow arguments:
        -profile                set to 'local' when running locally, set to 'singularity' when running on the HPC


        Mandatory pipeline arguments:
        --input_bed                     path to input file with fusionRNAs in bed format (0-based annotation)
        --index_bowtie          path to bowtie genome index directory
        --index_bowtie_name     the basename of the Bowtie index to be searched (the name of any of the index files up to but not including the final .1.ebwt / .rev.1.ebwt / ...)
        --index_fasta           path to directory that contains the fastahack genome and index file
        --index_fasta_name      the name of the fastahack genome file

        Optional pipeline arguments:
        --splice                        when set to 'yes' the input sequence will be spliced, when set to 'no' the input sequence will be unspliced
        --primer_settings       path to file with primer3 settings (see primer3 manual)
        --chrom_file            file containing all chromosome sizes (to validate bed file)
        --known_exons           bed file containing exon annotation
        --list_ENST                     file containing ENST numbers of canonical transcripts or transcripts of interest (this file can also be left empty)
        --primer3_diff          the minimum number of base pairs between the 3' ends of any two left primers (see also primer3 PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE)
        --primer3_nr            the number of primers designed by primer3; caution: setting this parameter to a large value will increase running time
        --min_tm                        minimum melt temperature of the primers (default: 58)
        --max_tm                        maximum melt temperature of the primers(default: 60)
        --opt_tm                        optimal melt temperature of the primers(default: 59)
        --diff_tm                       maximum difference in melt temperature between the primers(default: 2)
        --min_gc                        minimum GC contect of the primers (default: 30)
        --max_gc                        maximum GC contect of the primers(default: 80)
        --opt_gc                        optimal GC contect of the primers(default: 50)
        --amp_min                       minimum amplicon length (default: 60)
        --amp_max                       maximum amplicon length (default: 0)
        --temp_l                        the number of nucleotides on each side of the breakpoint that will be used for the template (example 150 => template of 300 nts in total)
        --spec_filter           when set to 'strict', a maximum of 2 MM is allowed; when set to 'loose', a maximum of 5 combined MM is allowed
        --snp_filter            when set to 'strict', no common SNPs are allowed in primer sequence; when set to 'loose', common SNPs are allowed in 5' first half of primer; when set to 'off', no filter is applied
        --snp_url                       when using a differente species than human, the correct SNP database url should be provided; alternatively, this paramater can be set to 'off' if no SNP database is available
        --upfront_filter        when set to 'yes', SNPs and secundary structures are avoided before primer design; when set to 'str', secundary structures are avoided before primer design; when set to 'snp', snp are avoided before primer design; when set to 'no', no filtering before primer design is performed
        --output_dir            path to directory where the output files will be saved
	
```


You can easily create your own profiles by modifying the nextflow.config file.

Nextflow keeps track of all the processes executed in your pipeline. If you want to restart the pipeline after a bug, add '-resume'. The execution of the processes that are not changed will be skipped and the cached result used instead.

Note: The pipeline results are cached by default in the directory $PWD/work. This folder can take of lot of disk space. If your are sure you won’t resume your pipeline execution, clean this folder periodically.

Note: If a fusionRNA is smaller than the requested template size, the template size is reduced to the fusionRNA size. Of note, if this 300-nucleotide template sequence includes an exon-intron boundary, the intronic region (which may not be part of the fusionRNA) is included. Some fusionRNAs effectively also include intronic sequences, and some BPs concatenate an exonic and an intronic sequence. 

### Running with nupack
If you have a copy of nupack (e.g. nupack-4.0.1.9) you can choose to run the pipeline with nupack instead of ViennaRNA. To do so you will have to build the docker image from the docker file included in /Docker. 
1. Uncomment the following lines in the docker file and make sure the version is correct.
```
#ADD ./nupack-4.0.1.9.zip /bin/
#RUN unzip nupack-4.0.1.9
#RUN python3 -m pip install -U nupack -f ./nupack-4.0.1.9/package
```
2. Build the image (in Docker folder)
```bash
docker build -t fusionprimerxl:both .
```
3. copy the include nextflow script (in bin) to the base folder
```
mv ./bin/FUSIONprimerXL_both.nf ./FUSIONprimerXL_both.nf
```
4. run the pipeline with the preferred program (ViennaRNA or Nupack)
```
nextflow run FUSIONprimerXL_both.nf -profile both --prediction_program Nupack
```

## 3. Step Running on the HPC (UGent)
<hr>

Nextflow version 20.10.0 is available on all clusters (swalot, skitty, victini, joltik, kirlia, doduo). The pipeline can be run through an interactive session. The pipeline can only run from the $VSC_SCRATCH_VO_USER directory.

```
qsub -I -l nodes=1:ppn=16 -l walltime=04:00:00
cd $VSC_SCRATCH_VO_USER/FUSIONprimerXL/
module load Nextflow/20.10.0
nextflow run FUSIONprimerXL.nf --help
```

## 4. Other species
<hr>

As default, FUSIONprimerXL designs primers for humans (GRCh38). To design primers for other species, the following files have to be provided and parsed through the corresponding parameters:

- A file containing the chromosome sizes (parameter **chrom_file**) (for example from: https://www.ncbi.nlm.nih.gov/grc/human/data)
- A [fastahack index](#Building-fastahack-index) (parameter **index_fasta**) and [Bowtie index](#Building-the-Bowtie-index)
- A SNP database link (parameter snp_url) (for example: http://hgdownload.soe.ucsc.edu/gbdb/hg38/snp/dbSnp155Common.bb)
- A file containing all ENST numbers of canonical transcripts or transcrits. This file can be generated by downloading the MANE file (https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/) and transforming it into a simple list using 00_C_generate_ENST_list.py (can be run in the Docker image).

#### Example references


| Species                | Reference       | Source                                                                                                                                                                                                                                                   |
| ---------------------- | --------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Xenopus tropicalis     | DNA             | [http://ftp.ensembl.org/pub/release-111/fasta/xenopus_tropicalis/dna/Xenopus_tropicalis.UCB_Xtro_10.0.dna.toplevel.fa.gz](http://ftp.ensembl.org/pub/release-111/fasta/xenopus_tropicalis/dna/Xenopus_tropicalis.UCB_Xtro_10.0.dna.toplevel.fa.gz)       |
|                        | cDNA            | [http://ftp.ensembl.org/pub/release-111/fasta/xenopus_tropicalis/cdna/Xenopus_tropicalis.UCB_Xtro_10.0.cdna.all.fa.gz](http://ftp.ensembl.org/pub/release-111/fasta/xenopus_tropicalis/cdna/Xenopus_tropicalis.UCB_Xtro_10.0.cdna.all.fa.gz)             |
|                        | ncRNA           | [http://ftp.ensembl.org/pub/release-111/fasta/xenopus_tropicalis/ncrna/Xenopus_tropicalis.UCB_Xtro_10.0.ncrna.fa.gz](http://ftp.ensembl.org/pub/release-111/fasta/xenopus_tropicalis/ncrna/Xenopus_tropicalis.UCB_Xtro_10.0.ncrna.fa.gz)                 |
|                        | exon annotation | [http://ftp.ensembl.org/pub/release-111/gtf/xenopus_tropicalis/Xenopus_tropicalis.UCB_Xtro_10.0.111.gtf.gz](http://ftp.ensembl.org/pub/release-111/gtf/xenopus_tropicalis/Xenopus_tropicalis.UCB_Xtro_10.0.111.gtf.gz)                                   |
| Mus musculus           | DNA             | [http://ftp.ensembl.org/pub/release-111/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz](http://ftp.ensembl.org/pub/release-111/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz)                             |
|                        | cDNA            | [http://ftp.ensembl.org/pub/release-111/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz](http://ftp.ensembl.org/pub/release-111/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz)                                                   |
|                        | ncRNA           | [http://ftp.ensembl.org/pub/release-111/fasta/mus_musculus/ncrna/Mus_musculus.GRCm39.ncrna.fa.gz](http://ftp.ensembl.org/pub/release-111/fasta/mus_musculus/ncrna/Mus_musculus.GRCm39.ncrna.fa.gz)                                                       |
|                        | exon annotation | [http://ftp.ensembl.org/pub/release-111/gtf/mus_musculus/Mus_musculus.GRCm39.111.gtf.gz](http://ftp.ensembl.org/pub/release-111/gtf/mus_musculus/Mus_musculus.GRCm39.111.gtf.gz)                                                                         |
| Caenorhabditis elegans | DNA             | [http://ftp.ensembl.org/pub/release-111/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz](http://ftp.ensembl.org/pub/release-111/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz) |
|                        | cDNA            | [http://ftp.ensembl.org/pub/release-111/fasta/caenorhabditis_elegans/cdna/Caenorhabditis_elegans.WBcel235.cdna.all.fa.gz](http://ftp.ensembl.org/pub/release-111/fasta/caenorhabditis_elegans/cdna/Caenorhabditis_elegans.WBcel235.cdna.all.fa.gz)       |
|                        | ncRNA           | [http://ftp.ensembl.org/pub/release-111/fasta/caenorhabditis_elegans/ncrna/Caenorhabditis_elegans.WBcel235.ncrna.fa.gz](http://ftp.ensembl.org/pub/release-111/fasta/caenorhabditis_elegans/ncrna/Caenorhabditis_elegans.WBcel235.ncrna.fa.gz)           |
|                        | exon annotation | [http://ftp.ensembl.org/pub/release-111/gtf/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.111.gtf.gz](http://ftp.ensembl.org/pub/release-111/gtf/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.111.gtf.gz)                             |
| Danio rerio            | DNA             | [http://ftp.ensembl.org/pub/release-111/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna.primary_assembly.fa.gz](http://ftp.ensembl.org/pub/release-111/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna.primary_assembly.fa.gz)                                 |
|                        | cDNA            | [http://ftp.ensembl.org/pub/release-111/fasta/danio_rerio/cdna/Danio_rerio.GRCz11.cdna.all.fa.gz](http://ftp.ensembl.org/pub/release-111/fasta/danio_rerio/cdna/Danio_rerio.GRCz11.cdna.all.fa.gz)                                                       |
|                        | ncRNA           | [http://ftp.ensembl.org/pub/release-111/fasta/danio_rerio/ncrna/Danio_rerio.GRCz11.ncrna.fa.gz](http://ftp.ensembl.org/pub/release-111/fasta/danio_rerio/ncrna/Danio_rerio.GRCz11.ncrna.fa.gz)                                                           |
|                        | exon annotation | [http://ftp.ensembl.org/pub/release-111/gtf/danio_rerio/Danio_rerio.GRCz11.111.gtf.gz](http://ftp.ensembl.org/pub/release-111/gtf/danio_rerio/Danio_rerio.GRCz11.111.gtf.gz)                                                                             |
| Rattus norvegicus      | DNA             | [http://ftp.ensembl.org/pub/release-111/fasta/rattus_norvegicus/dna/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa.gz](http://ftp.ensembl.org/pub/release-111/fasta/rattus_norvegicus/dna/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa.gz)                   |
|                        | cDNA            | [http://ftp.ensembl.org/pub/release-111/fasta/rattus_norvegicus/cdna/Rattus_norvegicus.mRatBN7.2.cdna.all.fa.gz](http://ftp.ensembl.org/pub/release-111/fasta/rattus_norvegicus/cdna/Rattus_norvegicus.mRatBN7.2.cdna.all.fa.gz)                         |
|                        | ncRNA           | [http://ftp.ensembl.org/pub/release-111/fasta/rattus_norvegicus/ncrna/Rattus_norvegicus.mRatBN7.2.ncrna.fa.gz](http://ftp.ensembl.org/pub/release-111/fasta/rattus_norvegicus/ncrna/Rattus_norvegicus.mRatBN7.2.ncrna.fa.gz)                             |
|                        | exon annotation | [http://ftp.ensembl.org/pub/release-111/gtf/rattus_norvegicus/Rattus_norvegicus.mRatBN7.2.111.gtf.gz](http://ftp.ensembl.org/pub/release-111/gtf/rattus_norvegicus/Rattus_norvegicus.mRatBN7.2.111.gtf.gz)                                               |




## 5. Nextflow tower
<hr>

[Nextflow tower](https://tower.nf/) can be used to monitor the pipeline while it's running.
```
nextflow run FUSIONprimerXL.nf -with-tower
```

When Nextflow tower is used in combination with the HPC, the nextflow version and tower access token should be indicated.
```
export NXF_VER="20.10.0"
export TOWER_ACCESS_TOKEN=your_token_here
```

## Cite

***PUBLICATION TO BE PLACED HERE*** 
