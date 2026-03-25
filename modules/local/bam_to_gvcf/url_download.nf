process URL_DOWNLOAD {
    tag "${sampleId}"

    conda "conda-forge::wget bioconda::samtools=1.17 bioconda::htslib=1.17"

    def cachedir = params.genotypes_cache ? file(params.genotypes_cache) : workDir
    storeDir cachedir / "bam_to_gvcf" / "bam"

    input:
    tuple val(sampleId), val(bam_url)

    output:
    tuple val(sampleId), path("${sampleId}.bam"), path("${sampleId}.bam.bai"), emit: bam

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
    elif wget -q --tries=3 --timeout=60 -O "!{sampleId}.bai" "${BAI_URL_2}" 2>/dev/null; then
        echo "INFO  Downloaded index from ${BAI_URL_2}"
    else
        echo "WARN  No remote BAI found — generating index with samtools"
        samtools index -@ !{task.cpus} "!{sampleId}.bam"
    fi

    # Normalize: rename .bai → .bam.bai if that is what was produced above
    # (covers the BAI_URL_2 branch and any future shell path that writes .bai)
    if [ -f "!{sampleId}.bai" ] && [ ! -f "!{sampleId}.bam.bai" ]; then
        mv "!{sampleId}.bai" "!{sampleId}.bam.bai"
    fi
    '''
}