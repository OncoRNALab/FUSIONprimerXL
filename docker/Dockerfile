# based on existing Docker image
FROM oncornalab/primerxl_circ:v0.27

# install python scripts
ADD scripts /usr/bin
RUN chmod +x /usr/bin/get_sec_str_temp_dmas.py /usr/bin/get_sec_str_amp_dmas.py /usr/bin/filter_dmas.py /usr/bin/split_primers_dmas.py /usr/bin/split_circRNAs_dmas.py /usr/bin/gather_output_dmas.py /usr/bin/get_SNPs_dmas.py /usr/bin/upfront_filter_dmas.py /usr/bin/validate_bed_dmas.py /usr/bin/summary_run_dmas.py /usr/bin/get_circ_seq_fastahack2_dmas.py /usr/bin/gtf_to_bed_dmas.py /usr/bin/generate_ENST_list_dmas.py
