process GENERATE_SAMPLESHEET {
    tag "samplesheet"

    publishDir "${params.outdir}/bam_to_gvcf/", mode: 'copy'

    input:
    val(gvcf_names)

    output:
    path("pgsc_calc_samplesheet.csv"), emit: csv

    script:
    def publish_dir = file(params.outdir).toAbsolutePath().resolve('bam_to_gvcf/gvcfs')
    def rows = gvcf_names.collect { name -> "cohort,${publish_dir}/${name},,," }.join('\n')
    """
    cat > pgsc_calc_samplesheet.csv <<'SAMPLESHEET'
sampleset,vcf_path,bfile_path,pfile_path,chrom
${rows}
SAMPLESHEET
    """
}