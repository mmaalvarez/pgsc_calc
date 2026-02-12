process FILTER_DBSNP {
    tag "filter-dbsnp"

    input:
    tuple path(cohort_ref), path(cohort_fai), path(cohort_dict)
    path(original_dbsnp)
    path(original_dbsnp_index)

    output:
    tuple path("dbsnp138_cohort_compatible.vcf.gz"),
          path("dbsnp138_cohort_compatible.vcf.gz.tbi"), emit: dbsnp

    shell:
    '''
    set -euo pipefail
    export TMPDIR=.

    bcftools view --threads !{task.cpus} \
      -R <(awk '{print $1"\\t0\\t"$2}' !{cohort_fai}) \
      !{original_dbsnp} \
      -Oz -o temp_filtered.vcf.gz

    bcftools reheader \
      -f !{cohort_fai} \
      temp_filtered.vcf.gz \
      -o dbsnp138_cohort_compatible.vcf.gz

    tabix -p vcf dbsnp138_cohort_compatible.vcf.gz
    rm -f temp_filtered.vcf.gz
    '''
}