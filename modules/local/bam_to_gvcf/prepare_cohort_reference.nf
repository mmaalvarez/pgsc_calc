process PREPARE_COHORT_REF {
    tag "cohort-reference"

    cachedir = params.genotypes_cache ? file(params.genotypes_cache) : workDir
    storeDir cachedir / "bam_to_gvcf" / "cohort_ref"
    
    input:
    tuple val(sampleId), path(one_bam), path(one_bai)
    path(base_reference)
    path(contig_map)

    output:
    tuple path('cohort_reference.fa'),
          path('cohort_reference.fa.fai'),
          path('cohort_reference.dict'), emit: reference

    shell:
    '''
    set -euo pipefail

    expand_reference.py \
        -b !{one_bam} \
        -f !{base_reference} \
        -m !{contig_map} \
        -o cohort_reference.fa

    samtools faidx cohort_reference.fa
    gatk --java-options "-Xmx3500m" CreateSequenceDictionary -R cohort_reference.fa -O cohort_reference.dict
    '''
}