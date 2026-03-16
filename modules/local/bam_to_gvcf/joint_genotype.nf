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
        -O cohort_raw_!{chromosome}.vcf.gz \
        -L !{chromosome} \
        --dbsnp !{dbsnp} \
        --include-non-variant-sites true

    ## ---- Step 3: Annotate with dbSNP and filter ----
        
    # 3.1 Pre-extract CHROM/POS/ALT for this chromosome only from dbSNP into a
    #     bgzipped TSV, then index it.  This avoids feeding the full dbSNP VCF
    #     to bcftools annotate, which caused cross-chromosome contamination
    
    bcftools query \\
        -r !{chromosome} \\
        -f '%CHROM\\t%POS\\t%ALT\\n' \\
        !{dbsnp} \\
        | bgzip -c > dbsnp_!{chromosome}.tsv.gz

    tabix -s 1 -b 2 -e 2 dbsnp_!{chromosome}.tsv.gz

    # 3.2 Write the extra INFO header line for the new DB_ALT field
    printf '##INFO=<ID=DB_ALT,Number=.,Type=String,Description="ALT allele(s) from dbSNP">\n' > db_alt.hdr

    # 3.3 Annotate with the chromosome-specific TSV.
    #     -c CHROM,POS,INFO/DB_ALT maps TSV columns 1/2/3 to those VCF fields.
    #     -r !{chromosome} is an extra safety guard to never emit other chroms.
    #     Every position in cohort_raw keeps its row; only positions that also
    #     exist in the TSV gain a DB_ALT=<allele> tag in their INFO column.
    
    bcftools annotate \\
        -a dbsnp_!{chromosome}.tsv.gz \\
        -c CHROM,POS,INFO/DB_ALT \\
        -h db_alt.hdr \\
        -r !{chromosome} \\
        -O z -o cohort_annotated_!{chromosome}.vcf.gz \\
        cohort_raw_!{chromosome}.vcf.gz

    # 3.4 awk pass over the annotated VCF:
    #   • header lines           → pass through unchanged
    #   • ALT is a real allele   → keep as-is (already a variant call)
    #   • ALT is "." or "*" AND DB_ALT present → replace ALT with dbSNP allele, keep
    #   • ALT is "." or "*" AND DB_ALT absent  → drop (hom-ref, not in dbSNP)
    #
    # Note: DB_ALT may be comma-separated for multi-allelic dbSNP entries; the
    # full string is placed in $5 so downstream tools can handle it normally.
    
    zcat cohort_annotated_!{chromosome}.vcf.gz | awk '
    BEGIN { FS="\\t"; OFS="\\t" }
    /^#/ { print; next }
    {
        if ($5 == "." || $5 == "*") {
            if ($8 ~ /DB_ALT=/) {
                # Ordered iteration avoids gawk unordered-hash surprises
                n = split($8, info, ";")
                db_alt = ""
                for (i = 1; i <= n; i++) {
                    if (info[i] ~ /^DB_ALT=/) {
                        sub(/^DB_ALT=/, "", info[i])
                        db_alt = info[i]
                        break
                    }
                }
                # && db_alt != $4 guard drops any site where the dbSNP ALT is identical to the REF — these are biologically meaningless entries in dbSNP
                if (db_alt != "" && db_alt != "." && db_alt != $4) {
                    $5 = db_alt
                    print
                }
                # else: DB_ALT tag present but empty/dot → still drop
            }
            # else: no DB_ALT tag → hom-ref not in dbSNP → drop
        } else {
            # Real ALT allele already called → keep unconditionally
            print
        }
    }' | bgzip -c > cohort_!{chromosome}.vcf.gz

    ## ---- Step 4: Index the joint-called VCF ----
    tabix -p vcf cohort_!{chromosome}.vcf.gz

    ## ---- Step 5: Clean up the intermediate combined gVCF ----
    rm -f cohort_combined.g.vcf.gz cohort_combined.g.vcf.gz.tbi \
        cohort_raw_!{chromosome}.vcf.gz cohort_raw_!{chromosome}.vcf.gz.tbi \
        cohort_annotated_!{chromosome}.vcf.gz db_alt.hdr
    '''
}