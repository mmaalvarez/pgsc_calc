process HAPLOTYPECALLER {
    tag "${sampleId}-${chromosome}-${scatter_idx}"

    conda "bioconda::gatk4=4.5.0.0 bioconda::bcftools=1.17 bioconda::htslib=1.17"

    input:
    tuple val(sampleId), path(normal_bam), path(normal_bai),
          path(reference), path(reference_fai), path(reference_dict),
          path(dbsnp), path(dbsnp_index),
          val(chromosome), val(scatter_idx), path(regions_interval)

    output:
    tuple val(sampleId), val(chromosome), val(scatter_idx),
          path("${sampleId}.${chromosome}.${scatter_idx}.g.vcf.gz"),
          path("${sampleId}.${chromosome}.${scatter_idx}.g.vcf.gz.tbi"), emit: gvcf

    shell:
    '''
    set -euo pipefail

    # 1. Call variants — raw gVCF for this sub-interval
    gatk --java-options "-Xmx!{Math.max((task.memory.toGiga() as int)-2,1)}G" HaplotypeCaller \
      -R !{reference} \
      -I !{normal_bam} \
      -O raw.g.vcf.gz \
      -L !{chromosome} \
      -L !{regions_interval} \
      --interval-set-rule INTERSECTION \
      --dbsnp !{dbsnp} \
      -ERC GVCF \
      --native-pair-hmm-threads !{task.cpus} \
      --create-output-variant-index true

    # 2. Filter: GQ>=20 and (MIN_)DP>=10, exclude REF="N"
    bcftools view \
      -i '(FORMAT/GQ>=20 & FORMAT/DP>=10 & TYPE!="ref") | (FORMAT/GQ>=20 & FORMAT/MIN_DP>=10 & TYPE="ref")' \
      raw.g.vcf.gz \
    | bcftools view -e 'REF="N"' -Oz \
      -o !{sampleId}.!{chromosome}.!{scatter_idx}.g.vcf.gz

    # 3. Index filtered gVCF
    tabix -p vcf !{sampleId}.!{chromosome}.!{scatter_idx}.g.vcf.gz

    # 4. Clean up
    rm -f raw.g.vcf.gz raw.g.vcf.gz.tbi
    '''
}