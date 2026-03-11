process HAPLOTYPECALLER {
    tag "${sampleId}-${chromosome}"

    input:
    tuple val(sampleId), path(normal_bam), path(normal_bai),
          path(reference), path(reference_fai), path(reference_dict),
          path(dbsnp), path(dbsnp_index),
          val(chromosome), path(regions_bed)

    output:
    tuple val(sampleId), val(chromosome),
          path("${sampleId}.${chromosome}.g.vcf.gz"),
          path("${sampleId}.${chromosome}.g.vcf.gz.tbi"), emit: gvcf

    shell:
    '''
    set -euo pipefail

    # 1. Call variants and produce raw gVCF
    gatk --java-options "-Xmx!{Math.max((task.memory.toGiga() as int)-2,1)}G" HaplotypeCaller \
      -R !{reference} \
      -I !{normal_bam} \
      -O !{sampleId}.!{chromosome}.raw.g.vcf.gz \
      -L !{chromosome} \
      -L !{regions_bed} \
      --interval-set-rule INTERSECTION \
      --dbsnp !{dbsnp} \
      -ERC GVCF \
      --native-pair-hmm-threads !{task.cpus} \
      --create-output-variant-index true

    # 2. Filter: keep sites with GQ>=20 and (MIN_)DP>=10, exclude REF="N"
    bcftools view -i '(FORMAT/GQ>=20 & FORMAT/DP>=10 & TYPE!="ref") | (FORMAT/GQ>=20 & FORMAT/MIN_DP>=10 & TYPE="ref")' !{sampleId}.!{chromosome}.raw.g.vcf.gz | \
      bcftools view -e 'REF="N"' -Oz -o !{sampleId}.!{chromosome}.g.vcf.gz

    # 3. Index filtered gVCF
    tabix -p vcf !{sampleId}.!{chromosome}.g.vcf.gz

    # 4. Clean up intermediate raw gVCF
    rm -f !{sampleId}.!{chromosome}.raw.g.vcf.gz \
          !{sampleId}.!{chromosome}.raw.g.vcf.gz.tbi
    '''
}