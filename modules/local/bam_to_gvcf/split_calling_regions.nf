process SPLIT_CALLING_REGIONS {
    tag "${chromosome}"

    input:
    tuple val(chromosome), path(callingRegionsGz)

    output:
    tuple val(chromosome), path("${chromosome}.bed"), emit: bed

    shell:
    '''
    set -euo pipefail
    gzip -dc !{callingRegionsGz} \
      | awk -v chr="!{chromosome}" 'BEGIN{OFS="\\t"} $1==chr {print}' \
      > !{chromosome}.bed
    '''
}