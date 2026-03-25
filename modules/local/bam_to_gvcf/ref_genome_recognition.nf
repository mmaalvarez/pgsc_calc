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

    # *** TEMPORARY DEBUG — remove once root cause is confirmed ***
    echo "[reco debug] sampleId      = '!{sampleId}'"      >&2
    echo "[reco debug] bamFile       = '!{bamFile}'"       >&2
    echo "[reco debug] baiFile       = '!{baiFile}'"       >&2
    echo "[reco debug] reference_dict= '!{reference_dict}'" >&2
    echo "[reco debug] work dir contents:" >&2
    ls -la >&2
    # *** END TEMPORARY DEBUG ***

    # BUG FIX 1: separate assignment from export so set -e catches failures
    needsRealign=$(detect_reference_genome.py -b '!{bamFile}' -d '!{reference_dict}')
    export needsRealign

    echo "[reco debug] needsRealign='${needsRealign}'" >&2
    '''
}