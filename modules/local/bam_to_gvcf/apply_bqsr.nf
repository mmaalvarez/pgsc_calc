process APPLY_BQSR {
    tag "${sampleId}_${chrom}"

    conda "bioconda::gatk4=4.5.0.0"

    input:
    tuple val(sampleId), path(dedup_bam), path(dedup_bai),
          path(reference), path(reference_fai), path(reference_dict),
          path(recal_table),
          val(chrom)

    output:
    tuple val(sampleId), val(chrom),
          path("${sampleId}_${chrom}_recal.bam"),
          path("${sampleId}_${chrom}_recal.bam.bai"), emit: bam

    shell:
    '''
    set -euo pipefail

    gatk --java-options "-Xmx!{Math.max((task.memory.toGiga() as int)-2,1)}G" ApplyBQSR \
      -R !{reference} \
      -I !{dedup_bam} \
      --bqsr-recal-file !{recal_table} \
      -L !{chrom} \
      -O !{sampleId}_!{chrom}_recal.bam

    gatk BuildBamIndex \
      -I !{sampleId}_!{chrom}_recal.bam \
      -O !{sampleId}_!{chrom}_recal.bam.bai
    '''
}