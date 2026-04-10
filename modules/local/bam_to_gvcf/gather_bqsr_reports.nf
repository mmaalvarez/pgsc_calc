process GATHER_BQSR_REPORTS {
    tag "${sampleId}"

    conda "bioconda::gatk4=4.5.0.0"

    input:
    tuple val(sampleId), path(tables)

    output:
    tuple val(sampleId),
          path("${sampleId}_merged_recal.table"), emit: table

    shell:
    '''
    set -euo pipefail

    input_args=""
    for tbl in !{tables instanceof List ? tables.join(' ') : tables}; do
        input_args="${input_args} -I ${tbl}"
    done

    gatk GatherBQSRReports \
      ${input_args} \
      -O !{sampleId}_merged_recal.table
    '''
}