process WGS_METRICS {
    tag "${sampleId}"

    publishDir "${params.outdir}/bam_to_gvcf/bam_metrics/", mode: 'copy'

    input:
    tuple val(sampleId), path(inputBam), path(inputBai),
          path(reference), path(reference_fai), path(reference_dict)

    output:
    tuple val(sampleId), path("${sampleId}_wgsMetrics.txt"), emit: metrics

    shell:
    '''
    set -euo pipefail
    gridssJarPath=$(locate_gridss_jar.sh)

    java -Xmx!{Math.max((task.memory.toGiga() as int)-1,1)}G \
      -Dsamjdk.use_async_io_read_samtools=true \
      -Dsamjdk.use_async_io_write_samtools=true \
      -Dsamjdk.use_async_io_write_tribble=true \
      -Dsamjdk.buffer_size=4194304 \
      -cp $gridssJarPath \
      picard.cmdline.PicardCommandLine CollectWgsMetrics \
      REFERENCE_SEQUENCE=!{reference} \
      INPUT=!{inputBam} \
      OUTPUT=!{sampleId}_wgsMetrics.txt \
      MINIMUM_MAPPING_QUALITY=20 \
      MINIMUM_BASE_QUALITY=10 \
      COVERAGE_CAP=250
    '''
}