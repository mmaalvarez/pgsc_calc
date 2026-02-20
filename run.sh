#!/bin/bash

## load modules, and set root dir names and partition names, depending on cluster ('executor')

hostname=`echo $HOSTNAME | cut -d"." -f1 | cut -d"-" -f1`

if [[ "$hostname" == "" ]]; then

    ## LUCIA cluster

elif [[ "$hostname" == "lsflogin" || "$hostname" == "hpc" ]]; then

    ## GEL cluster

    module load nextflow/24.04.2-with-plugins
    module load singularity/4.1.1
    module load graphviz/12.1.1

    export root_dir="/nas/weka.gel.zone"
    export work_dir=""$root_dir"/re_gecip/cancer_pan/malvarez/lucia/"
    export executor="lsf"
    export partition_fast_short="short"
    export partition_slow_long="medium"
    export partition_slowest_unlimited="long"
    export project_group="-P re_gecip_cancer_pan"
    export TMPDIR="/re_scratch/re_gecip/cancer_pan/malvarez"

elif [[ "$hostname" == "fsupeksvr" ]]; then

    ## AGENDAS cluster

    conda activate nextflow

    export root_dir="/g"
    export work_dir=""$root_dir"/strcombio/fsupek_data/users/malvarez/projects/lucia/"
    ref_dir=""$root_dir"/strcombio/fsupek_fisher/malvarez/GWSPipeline/utils/"
    export executor="slurm"
    export partition_fast_short="normal_prio"
    export partition_slow_long="normal_prio_long"
    export partition_slowest_unlimited="normal_prio_unlim"
    export project_group=""

else
    echo "ERROR: HOSTNAME is not known: '`echo $HOSTNAME | cut -d"." -f1`'"
fi


# if the input is a BAM files (paths or download IDs) table (TSV with a 'bamFile' column), run:
nextflow run mmaalvarez/pgsc_calc \
    --input /path/to/BAM_paths_or_ids.tsv \
    --pgs_id PGS000740 \
    --target_build GRCh37 \
    --bam2gvcf_fasta "$ref_dir"/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
    --bam2gvcf_dict "$ref_dir"/GCA_000001405.15_GRCh38_no_alt_analysis_set.dict \
    --bam2gvcf_dbsnp "$ref_dir"/Homo_sapiens_assembly38.dbsnp138.vcf.gz \
    --bam2gvcf_calling_regions "$ref_dir"/GRCh38.d1.vd1.fa.bed.gz \
    --bam2gvcf_contig_map "$ref_dir"/contig_mappings.txt \
    --gdc_token "$ref_dir"/gdc_token.txt \
    -profile singularity \
    -resume

# if the input is already a gVCF paths table (CSV with the columns required by pgsc_calc) it only needs:
#nextflow run mmaalvarez/pgsc_calc --input /path/to/gVCF_paths.csv --pgs_id PGS000740 --target_build GRCh37 -profile singularity -resume
