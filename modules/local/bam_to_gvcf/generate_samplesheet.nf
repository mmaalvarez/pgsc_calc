process GENERATE_SAMPLESHEET {
    tag "samplesheet"

    publishDir "${params.outdir}/bam_to_gvcf/", mode: 'copy'

    input:
    val(vcf_name)

    output:
    path("pgsc_calc_samplesheet.csv"), emit: csv

    script:
    // point to the published location of the joint-called VCF
    def publish_dir = file(params.outdir).toAbsolutePath().resolve('bam_to_gvcf/joint_called')
    """
    cat > pgsc_calc_samplesheet.csv <<'SAMPLESHEET'
sampleset,vcf_path,bfile_path,pfile_path,chrom
cohort,${publish_dir}/${vcf_name},,,
SAMPLESHEET
    """
}