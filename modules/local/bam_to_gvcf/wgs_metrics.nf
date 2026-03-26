process WGS_METRICS {
    tag "${sampleId}"

    conda "bioconda::gridss=2.11.1"

    publishDir "${params.outdir}/bam_to_gvcf/bam_metrics/", mode: 'copy'

    input:
    tuple val(sampleId), path(inputBam), path(inputBai),
          path(reference), path(reference_fai), path(reference_dict)

    output:
    tuple val(sampleId), path("${sampleId}_wgsMetrics.txt"), emit: metrics

    shell:
    // *** CHANGED: replaced fragile `locate_gridss_jar.sh` + manual `java -cp` invocation
    // with `gatk CollectWgsMetrics`, which is already available from bioconda::gatk4=4.5.0.0
    // (used by DEDUP_BQSR, HAPLOTYPECALLER, JOINT_GENOTYPE, PREPARE_COHORT_REF). ***
    '''
    set -euo pipefail
    gatk --java-options "-Xmx!{Math.max((task.memory.toGiga() as int)-1,1)}G" \
      CollectWgsMetrics \
      REFERENCE_SEQUENCE=!{reference} \
      INPUT=!{inputBam} \
      OUTPUT=!{sampleId}_wgsMetrics.txt \
      MINIMUM_MAPPING_QUALITY=20 \
      MINIMUM_BASE_QUALITY=10 \
      COVERAGE_CAP=250
    '''
}