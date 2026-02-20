process JOINT_GENOTYPE {
    tag "joint_genotype"

    publishDir "${params.outdir}/bam_to_gvcf/joint_called/", mode: 'copy'

    input:
    path gvcfs
    path tbis
    tuple path(ref_fasta), path(ref_fai), path(ref_dict)
    tuple path(dbsnp), path(dbsnp_tbi)

    output:
    tuple path("cohort_joint_called.vcf.gz"),
          path("cohort_joint_called.vcf.gz.tbi"), emit: vcf

    shell:
    '''
    set -euo pipefail

    ## ---- build --variant arguments from the staged gVCFs ----
    variant_args=""
    for gvcf in !{gvcfs instanceof List ? gvcfs.join(' ') : gvcfs}; do
        variant_args="${variant_args} --variant ${gvcf}"
    done

    ## ---- Step 1: CombineGVCFs → multi-sample gVCF ----
    gatk --java-options "-Xmx!{Math.max((task.memory.toGiga() as int) - 4, 2)}G" CombineGVCFs \
        -R !{ref_fasta} \
        ${variant_args} \
        -O cohort_combined.g.vcf.gz \
        --dbsnp !{dbsnp}

    ## ---- Step 2: GenotypeGVCFs (include non-variant sites) ----
    gatk --java-options "-Xmx!{Math.max((task.memory.toGiga() as int) - 4, 2)}G" GenotypeGVCFs \
        -R !{ref_fasta} \
        -V cohort_combined.g.vcf.gz \
        -O cohort_joint_called.vcf.gz \
        --dbsnp !{dbsnp} \
        --include-non-variant-sites true

    ## ---- Index the joint-called VCF ----
    gatk --java-options "-Xmx2G" IndexFeatureFile \
        -I cohort_joint_called.vcf.gz

    ## ---- Clean up the intermediate combined gVCF ----
    rm -f cohort_combined.g.vcf.gz cohort_combined.g.vcf.gz.tbi
    '''
}