process MERGE_GVCFS {
    tag "${sampleId}_${chrom}"

    conda "bioconda::gatk4=4.5.0.0 bioconda::htslib=1.17"

    input:
    tuple val(sampleId), val(chrom), path(gvcfs), path(tbis)

    output:
    tuple val(sampleId), val(chrom),
          path("${sampleId}.${chrom}.g.vcf.gz"),
          path("${sampleId}.${chrom}.g.vcf.gz.tbi"), emit: gvcf

    shell:
    '''
    set -euo pipefail

    input_args=""
    for g in !{gvcfs instanceof List ? gvcfs.join(' ') : gvcfs}; do
        input_args="${input_args} -I ${g}"
    done

    gatk MergeVcfs \
      ${input_args} \
      -O !{sampleId}.!{chrom}.g.vcf.gz
    '''
}