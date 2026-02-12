process REF_GENOME_RECOGNITION {
    tag "${sampleId}"

    input:
    tuple val(sampleId), path(bamFile), path(baiFile)
    path(reference_dict)

    output:
    tuple val(sampleId), env(isGRCh38), path(bamFile), path(baiFile), emit: result

    shell:
    '''
    set -euo pipefail
    export isGRCh38=$(detect_reference_genome.py -b !{bamFile} -d !{reference_dict})
    '''
}