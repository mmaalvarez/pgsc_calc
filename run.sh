#!/bin/bash

mkdir -p log/

## if the input is already a multisample-VCF path table (CSV with the columns required by pgsc_calc) it only needs:
nextflow -log $PWD/log/nextflow.log run mmaalvarez/pgsc_calc \
    --input /path/to/input_table.csv \
    --pgs_id PGS000740 \
    --target_build GRCh37 \
    -profile singularity \
    -resume #\

## if the input is a single-sample-gVCF paths table (TSV containing a column with the 'gvcfFile' header), add:
#    --bam2gvcf_fasta /path/to/refDir/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
#    --bam2gvcf_dict /path/to/refDir/GCA_000001405.15_GRCh38_no_alt_analysis_set.dict \
#    --bam2gvcf_dbsnp /path/to/refDir/Homo_sapiens_assembly38.dbsnp138.vcf.gz #\

## and if the input is a BAM files (paths or download IDs) table (TSV containing a column with the 'bamFile' header), add also:
#    --bam2gvcf_calling_regions /path/to/refDir/GRCh38.d1.vd1.fa.bed.gz \
#    --bam2gvcf_contig_map /path/to/refDir/contig_mappings.txt \
#    --gdc_token /path/to/refDir/gdc_token.txt
