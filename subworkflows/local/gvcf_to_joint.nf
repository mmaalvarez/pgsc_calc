include { JOINT_GENOTYPE       } from '../../modules/local/bam_to_gvcf/joint_genotype'
include { GENERATE_SAMPLESHEET } from '../../modules/local/bam_to_gvcf/generate_samplesheet'

workflow GVCF_TO_JOINT {

    take:
    ch_sample_table   // channel: path to TSV with gvcfFile column

    main:

    // ---- Constants --------------------------------------------------------
    def CHROMS = [
        'chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10',
        'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20',
        'chr21','chr22'
    ]

    // ---- Reference files from params --------------------------------------
    ch_ref_tuple = Channel.value([
        file(params.bam2gvcf_fasta,              checkIfExists: true),
        file("${params.bam2gvcf_fasta}.fai",     checkIfExists: true),
        file(params.bam2gvcf_dict,               checkIfExists: true)
    ])

    ch_dbsnp_tuple = Channel.value([
        file(params.bam2gvcf_dbsnp,              checkIfExists: true),
        file("${params.bam2gvcf_dbsnp}.tbi",     checkIfExists: true)
    ])

    // ---- Parse sample table -----------------------------------------------
    //  Expects a TAB-separated file with (at least) a  gvcfFile  column.
    //  Each gVCF must have a .tbi sitting alongside.
    ch_gvcf_pairs = ch_sample_table
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            if (!row.gvcfFile)
                error "Missing 'gvcfFile' column in row: ${row}"
            def gvcf = file(row.gvcfFile as String, checkIfExists: true)
            def tbi  = file("${row.gvcfFile}.tbi",  checkIfExists: true)
            tuple(gvcf, tbi)
        }

    // Collect all whole-genome gVCFs and TBIs, wrapping each collected list
    // inside another list so that combine preserves them as nested lists
    // (rather than flattening into the outer tuple).
    ch_gvcf_files = ch_gvcf_pairs.map { gvcf, tbi -> gvcf }.collect().map { [it] }
    ch_tbi_files  = ch_gvcf_pairs.map { gvcf, tbi -> tbi  }.collect().map { [it] }

    // ---- Per-chromosome joint genotyping ----------------------------------
    //  Every chromosome gets ALL samples' gVCFs; the -L flag inside
    //  JOINT_GENOTYPE restricts CombineGVCFs + GenotypeGVCFs to that chrom.
    Channel.of(*CHROMS)
        .combine(ch_gvcf_files)
        .combine(ch_tbi_files)
        .combine(ch_ref_tuple)
        .combine(ch_dbsnp_tuple)
        .set { ch_jg_input }
    // ch_jg_input: [chrom, [gvcfs], [tbis], ref, fai, dict, dbsnp, dbsnp_tbi]

    JOINT_GENOTYPE(ch_jg_input)

    // ---- Generate pgsc_calc samplesheet (one row per chromosome) ----------
    JOINT_GENOTYPE.out.vcf
        .map { chrom, vcf, tbi -> vcf }
        .collect()
        .set { ch_all_vcfs }

    GENERATE_SAMPLESHEET(ch_all_vcfs)

    emit:
    samplesheet = GENERATE_SAMPLESHEET.out.csv   // path to multi-row CSV
    joint_vcfs  = JOINT_GENOTYPE.out.vcf          // tuple(chrom, vcf, tbi) per chromosome
}