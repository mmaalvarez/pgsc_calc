process GATHER_GVCFS {
    tag "${sampleId}"

    publishDir "${params.outdir}/bam_to_gvcf/gvcfs/", mode: 'copy'

    input:
    tuple val(sampleId), val(chromosome_list), path(gvcf_list), path(tbi_list)

    output:
    tuple val(sampleId),
          path("${sampleId}.g.vcf.gz"),
          path("${sampleId}.g.vcf.gz.tbi"), emit: gvcf

    shell:
    '''
    set -euo pipefail

    input_args=""
    for chrom in !{chromosome_list.join(' ')}; do
      input_args="${input_args} -I !{sampleId}.${chrom}.g.vcf.gz"
    done

    gatk --java-options "-Xmx!{Math.max((task.memory.toGiga() as int)-2,1)}G" GatherVcfs \
      ${input_args} \
      -O !{sampleId}.g.vcf.gz \
      --CREATE_INDEX true \
      --REORDER_INPUT_BY_FIRST_VARIANT true
    '''
}