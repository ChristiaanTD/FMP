#!/bin/bash

#Due to temp-memory limitations analysis were not run on all samples simultaneously. This script for example only includes sample 23 & 24 but of course is applicable to all samples upon tweaking the input dir parameter!

cd "/data/main_storage/"

nextflow run nf-core/rnaseq -r 3.14.0 \
-resume \
-profile docker \
--input samplesheet_23_24.csv \
--outdir './results_23_24' \
--reads  './raw_files/incl_umi/*{_R1,_R3}_with_UMIs.fastq.gz' \
--fasta './essentials/Homo_sapiens.GRCh38.dna.primary_assembly.fa' \
--gtf  './essentials/Homo_sapiens.GRCh38.112.gtf' \
--with_umi \
--skip_umi_extract \
--umitools_umi_separator ":" \
--umitools_grouping_method percentile \
--trimmer fastp \
--igenomes_ignore \
--email t.d.dormans@amsterdamumc.nl \
-c nextflow.config

