process GATHER_BAM_FILES {
    tag "${sampleId}"

    conda "bioconda::gatk4=4.5.0.0"

    // Preserve the original storeDir caching for the final BAM
    def cachedir = params.genotypes_cache ? file(params.genotypes_cache) : workDir
    storeDir cachedir / "bam_to_gvcf" / "dedup_bqsr"

    input:
    tuple val(sampleId), path(bams), path(bais)

    output:
    tuple val(sampleId),
          path("${sampleId}_final.bam"),
          path("${sampleId}_final.bam.bai"), emit: bam

    shell:
    '''
    set -euo pipefail

    ## Build -I arguments in reference chromosome order (chr1 … chr22)
    input_args=""
    for i in $(seq 1 22); do
        f="!{sampleId}_chr${i}_recal.bam"
        if [ -f "${f}" ]; then
            input_args="${input_args} -I ${f}"
        fi
    done

    gatk GatherBamFiles \
      ${input_args} \
      -O !{sampleId}_final.bam

    gatk BuildBamIndex \
      -I !{sampleId}_final.bam \
      -O !{sampleId}_final.bam.bai
    '''
}