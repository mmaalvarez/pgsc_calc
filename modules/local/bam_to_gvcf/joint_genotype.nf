process JOINT_GENOTYPE {
    tag "joint_genotype_${chromosome}"

    publishDir "${params.outdir}/bam_to_gvcf/joint_called/", mode: 'copy'

    input:
    tuple val(chromosome), path(gvcfs), path(tbis),
          path(ref_fasta), path(ref_fai), path(ref_dict),
          path(dbsnp), path(dbsnp_tbi)

    output:
    tuple val(chromosome),
          path("cohort_${chromosome}.vcf.gz"),
          path("cohort_${chromosome}.vcf.gz.tbi"), emit: vcf

    shell:
    '''
    set -euo pipefail

    ## ---- build --variant arguments from the staged gVCFs ----
    variant_args=""
    for gvcf in !{gvcfs instanceof List ? gvcfs.join(' ') : gvcfs}; do
        variant_args="${variant_args} --variant ${gvcf}"
    done

    ## ---- Step 1: CombineGVCFs → multi-sample gVCF for this chromosome ----
    gatk --java-options "-Xmx!{Math.max((task.memory.toGiga() as int) - 4, 2)}G" CombineGVCFs \
        -R !{ref_fasta} \
        ${variant_args} \
        -L !{chromosome} \
        -O cohort_combined.g.vcf.gz \
        --dbsnp !{dbsnp}

    ## ---- Step 2: GenotypeGVCFs (restricted to this chromosome) ----
    gatk --java-options "-Xmx!{Math.max((task.memory.toGiga() as int) - 4, 2)}G" GenotypeGVCFs \
        -R !{ref_fasta} \
        -V cohort_combined.g.vcf.gz \
        -O cohort_!{chromosome}.vcf.gz \
        -L !{chromosome} \
        --dbsnp !{dbsnp} \
        --include-non-variant-sites true

    ## ---- Index the joint-called VCF ----
    gatk --java-options "-Xmx2G" IndexFeatureFile \
        -I cohort_!{chromosome}.vcf.gz

    ## ---- Clean up the intermediate combined gVCF ----
    rm -f cohort_combined.g.vcf.gz cohort_combined.g.vcf.gz.tbi
    '''
}