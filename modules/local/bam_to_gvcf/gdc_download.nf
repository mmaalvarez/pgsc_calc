process GDC_DOWNLOAD {
    tag "${gdc_id}"
    
    cachedir = params.genotypes_cache ? file(params.genotypes_cache) : workDir
    storeDir cachedir / "bam_to_gvcf" / "bam"

    input:
    val(gdc_id)

    output:
    tuple val(gdc_id), path("*.bam"), path("*.bai"), emit: bam

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
    find ./step_output_dir -name "*.bam" -exec mv -t . {} +
    find ./step_output_dir -name "*.bai" -exec mv -t . {} +
    '''
}