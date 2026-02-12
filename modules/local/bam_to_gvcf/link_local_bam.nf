process LINK_LOCAL_BAM {
    tag "${sampleId}"

    input:
    tuple val(sampleId), path(bam_file), path(bai_file)

    output:
    tuple val(sampleId), path("*.bam"), path("*.bai"), emit: bam

    script:
    """
    set -euo pipefail
    ls *.bam *.bai
    """
}