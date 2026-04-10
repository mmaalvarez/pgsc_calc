process BASE_RECALIBRATOR {
    tag "${sampleId}_${chrom}"

    conda "bioconda::gatk4=4.5.0.0"

    input:
    tuple val(sampleId), path(dedup_bam), path(dedup_bai),
          path(reference), path(reference_fai), path(reference_dict),
          path(dbsnp), path(dbsnp_index),
          val(chrom)

    output:
    tuple val(sampleId),
          path("${sampleId}_${chrom}_recal.table"), emit: table

    shell:
    '''
    set -euo pipefail

    gatk --java-options "-Xmx!{Math.max((task.memory.toGiga() as int)-2,1)}G" BaseRecalibrator \
      -R !{reference} \
      -I !{dedup_bam} \
      --known-sites !{dbsnp} \
      -L !{chrom} \
      -O !{sampleId}_!{chrom}_recal.table
    '''
}