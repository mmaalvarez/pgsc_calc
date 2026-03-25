process GDC_DOWNLOAD {
    tag "${sampleId}"                                   // *** CHANGED: was gdc_id ***

    conda "bioconda::gdc-client=1.6.1 bioconda::samtools=1.17"  // *** CHANGED: added samtools ***

    def cachedir = params.genotypes_cache ? file(params.genotypes_cache) : workDir  // *** CHANGED: added def ***
    storeDir cachedir / "bam_to_gvcf" / "bam"

    input:
    tuple val(sampleId), val(gdc_id)                   // *** CHANGED: now receives (sampleId, uuid) tuple ***

    output:
    // *** CHANGED: explicit names instead of globs; globs return lists, not single paths ***
    tuple val(sampleId), path("${sampleId}.bam"), path("${sampleId}.bam.bai"), emit: bam

    shell:
    '''
    set -euo pipefail
    mkdir -p step_output_dir
    gdc-client download \
        -n !{task.cpus} \
        -t !{params.gdc_token} \
        -d ./step_output_dir \
        --retry-amount 10 \
        !{gdc_id}

    # GDC places the download in a UUID-named subdirectory; locate and rename the BAM
    downloaded_bam=$(find ./step_output_dir -name "*.bam" | head -1)
    if [ -z "$downloaded_bam" ]; then
        echo "ERROR: gdc-client produced no BAM under step_output_dir" >&2
        exit 1
    fi
    mv "$downloaded_bam" "!{sampleId}.bam"

    # Prefer .bam.bai; fall back to .bai; generate with samtools if absent
    downloaded_bai=$(find ./step_output_dir -name "*.bam.bai" | head -1)
    if [ -z "$downloaded_bai" ]; then
        downloaded_bai=$(find ./step_output_dir -name "*.bai" | head -1)
    fi

    if [ -n "$downloaded_bai" ]; then
        mv "$downloaded_bai" "!{sampleId}.bam.bai"
    else
        echo "WARN  No BAI found in GDC download — generating with samtools" >&2
        samtools index -@ !{task.cpus} "!{sampleId}.bam"
    fi
    '''
}