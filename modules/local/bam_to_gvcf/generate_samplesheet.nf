process GENERATE_SAMPLESHEET {
    tag "samplesheet"

    publishDir "${params.outdir}/bam_to_gvcf/", mode: 'copy'

    input:
    path(vcf_files)

    output:
    path("pgsc_calc_samplesheet.csv"), emit: csv

    shell:
    '''
    set -euo pipefail

    echo "sampleset,path_prefix,chrom,format" > pgsc_calc_samplesheet.csv

    for vcf in cohort_*.vcf.gz; do
        [ -f "$vcf" ] || continue

        real_vcf=$(readlink -f "$vcf")

        # Extract chromosome from filename: cohort_chr1.vcf.gz -> chr1
        chrom=$(basename "$vcf" .vcf.gz | sed 's/^cohort_//')

        # Strip .vcf.gz to get path_prefix (SamplesheetParser re-appends the extension)
        prefix="${real_vcf%.vcf.gz}"

        echo "cohort,${prefix},${chrom},vcf" >> pgsc_calc_samplesheet.csv
    done

    # Sort rows by chromosome (natural/version sort) for reproducibility
    {
        head -1 pgsc_calc_samplesheet.csv
        tail -n +2 pgsc_calc_samplesheet.csv | sort -t',' -k3 -V
    } > tmp_sorted.csv
    mv tmp_sorted.csv pgsc_calc_samplesheet.csv
    '''
}