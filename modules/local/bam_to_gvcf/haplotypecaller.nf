process HAPLOTYPECALLER {
    tag "${sampleId}-${chromosome}"

    input:
    tuple val(sampleId), path(normal_bam), path(normal_bai),
          path(reference), path(reference_fai), path(reference_dict),
          path(dbsnp), path(dbsnp_index),
          val(chromosome), path(regions_bed)

    output:
    tuple val(sampleId), val(chromosome),
          path("${sampleId}.${chromosome}.g.vcf.gz"),
          path("${sampleId}.${chromosome}.g.vcf.gz.tbi"), emit: gvcf

    shell:
    '''
    set -euo pipefail

    gatk --java-options "-Xmx!{Math.max((task.memory.toGiga() as int)-2,1)}G" HaplotypeCaller \
      -R !{reference} \
      -I !{normal_bam} \
      -O !{sampleId}.!{chromosome}.g.vcf.gz \
      -L !{chromosome} \
      -L !{regions_bed} \
      --interval-set-rule INTERSECTION \
      --dbsnp !{dbsnp} \
      -ERC GVCF \
      --native-pair-hmm-threads !{task.cpus} \
      --create-output-variant-index true
    '''
}