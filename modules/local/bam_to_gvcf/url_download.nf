process URL_DOWNLOAD {
    tag "${sampleId}"

    input:
    tuple val(sampleId), val(bam_url)

    output:
    tuple val(sampleId), path("*.bam"), path("*.bai"), emit: bam

    shell:
    '''
    set -euo pipefail

    # Download BAM
    wget -q --tries=5 --timeout=120 -O "!{sampleId}.bam" "!{bam_url}"

    # Try the two common BAI naming conventions; fall back to samtools index
    BAI_URL_1="!{bam_url}.bai"                                         # foo.bam.bai
    BAI_URL_2="$(echo '!{bam_url}' | sed 's/\\.bam$/.bai/')"          # foo.bai

    if wget -q --tries=3 --timeout=60 -O "!{sampleId}.bam.bai" "${BAI_URL_1}" 2>/dev/null; then
        echo "INFO  Downloaded index from ${BAI_URL_1}"
    elif wget -q --tries=3 --timeout=60 -O "!{sampleId}.bam.bai" "${BAI_URL_2}" 2>/dev/null; then
        echo "INFO  Downloaded index from ${BAI_URL_2}"
    else
        echo "WARN  No remote BAI found — generating index with samtools"
        samtools index -@ !{task.cpus} "!{sampleId}.bam"
    fi
    '''
}
