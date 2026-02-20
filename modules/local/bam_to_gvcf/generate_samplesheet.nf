process GENERATE_SAMPLESHEET {
    tag "samplesheet"

    publishDir "${params.outdir}/bam_to_gvcf/", mode: 'copy'

    input:
    tuple path(vcf_file), path(vcf_tbi)

    output:
    path("pgsc_calc_samplesheet.csv"), emit: csv

    shell:
    '''
    set -euo pipefail

    # Resolve the staged symlink to the real VCF path in the work directory
    real_vcf=$(readlink -f !{vcf_file})

    # SamplesheetParser.getFilePaths() appends .vcf / .vcf.gz / .vcf.zst,
    # so path_prefix must be the path WITHOUT the file extension
    prefix="$real_vcf"
    if [[ "$real_vcf" == *.g.vcf.gz ]]; then
        prefix="${real_vcf%.g.vcf.gz}"
    elif [[ "$real_vcf" == *.g.vcf.zst ]]; then
        prefix="${real_vcf%.g.vcf.zst}"
    elif [[ "$real_vcf" == *.vcf.gz ]]; then
        prefix="${real_vcf%.vcf.gz}"
    elif [[ "$real_vcf" == *.vcf.zst ]]; then
        prefix="${real_vcf%.vcf.zst}"
    elif [[ "$real_vcf" == *.vcf ]]; then
        prefix="${real_vcf%.vcf}"
    fi

    cat > pgsc_calc_samplesheet.csv <<SAMPLESHEET
sampleset,path_prefix,chrom,format
cohort,${prefix},,vcf
SAMPLESHEET
    '''
}