process REF_GENOME_RECOGNITION {
    tag "${sampleId}"

    conda "bioconda::pysam=0.22.0"

    input:
    tuple val(sampleId), path(bamFile), path(baiFile)
    path(reference_dict)

    output:
    tuple val(sampleId), env(needsRealign), path("${sampleId}.bam"), path("${sampleId}.bam.bai"), emit: result

    shell:
    '''
    set -euo pipefail

    target_bam="!{sampleId}.bam"
    target_bai="${target_bam}.bai"

    # *** CHANGED: was missing BAM normalization entirely, and BAI normalization only
    # handled the special case "sampleId.bai" while ignoring the staged !{baiFile}.
    # Fix: create sampleId-based symlinks for both BAM and BAI whenever the staged
    # filename differs from the expected target name. This covers every upstream
    # source (LINK_LOCAL_BAM with its original filename, URL_DOWNLOAD, GDC_DOWNLOAD). ***

    # Normalise BAM filename
    [ -f "${target_bam}" ] || ln -sf "!{bamFile}" "${target_bam}"

    # Normalise BAI filename (handles .bam.bai, .bai, or any other staged name)
    [ -f "${target_bai}" ] || ln -sf "!{baiFile}" "${target_bai}"

    # Assign before export so set -e can catch a non-zero exit from the script
    needsRealign=$(detect_reference_genome.py -b "${target_bam}" -d '!{reference_dict}')
    export needsRealign
    '''
}