process FLAGSTAT {
    tag "${sampleId}"

    publishDir "${params.outdir}/bam_to_gvcf/bam_metrics/", mode: 'copy'

    input:
    tuple val(sampleId), path(inputBam), path(inputBai)

    output:
    tuple val(sampleId), path("${sampleId}.flagstat"), emit: stats

    shell:
    '''
    set -euo pipefail
    sambamba flagstat -t !{task.cpus} !{inputBam} > !{sampleId}.flagstat
    '''
}