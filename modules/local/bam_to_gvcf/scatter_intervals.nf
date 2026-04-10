process SCATTER_INTERVALS {
    tag "${chrom}"

    conda "bioconda::gatk4=4.5.0.0"

    input:
    tuple val(chrom), path(bed)
    path(ref_fasta)
    path(ref_fai)
    path(ref_dict)
    val(scatter_count)

    output:
    tuple val(chrom), path("scattered/*"), emit: intervals

    shell:
    '''
    set -euo pipefail
    mkdir -p scattered

    gatk SplitIntervals \
        -R !{ref_fasta} \
        -L !{bed} \
        --scatter-count !{scatter_count} \
        --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION \
        -O scattered
    '''
}