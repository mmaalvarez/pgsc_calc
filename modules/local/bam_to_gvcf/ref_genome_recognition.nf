process REF_GENOME_RECOGNITION {
    tag "${sampleId}"

    conda "bioconda::pysam=0.22.0"

    input:
    tuple val(sampleId), path(bamFile), path(baiFile)
    path(reference_dict)

    output:
    tuple val(sampleId), env(needsRealign), path(bamFile), path(baiFile), emit: result

    shell:
    '''
    set -euo pipefail
    export needsRealign=$(detect_reference_genome.py -b !{bamFile} -d !{reference_dict})
    '''
}