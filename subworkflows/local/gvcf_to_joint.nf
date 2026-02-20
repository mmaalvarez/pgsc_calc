include { JOINT_GENOTYPE       } from '../../modules/local/bam_to_gvcf/joint_genotype'
include { GENERATE_SAMPLESHEET } from '../../modules/local/bam_to_gvcf/generate_samplesheet'

workflow GVCF_TO_JOINT {

    take:
    ch_sample_table   // channel: path to TSV with gvcfFile column

    main:

    // ---- Reference files from params --------------------------------------
    //  Only the FASTA (+.fai), dict, and dbSNP (+.tbi) are needed;
    //  BWA indices, calling regions, etc. are NOT required for this path.
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

    ch_gvcf_files = ch_gvcf_pairs.map { gvcf, tbi -> gvcf }.collect()
    ch_tbi_files  = ch_gvcf_pairs.map { gvcf, tbi -> tbi  }.collect()

    // ---- Joint genotype: CombineGVCFs + GenotypeGVCFs ---------------------
    JOINT_GENOTYPE(
        ch_gvcf_files,
        ch_tbi_files,
        ch_ref_tuple,
        ch_dbsnp_tuple
    )

    // ---- Generate pgsc_calc-compatible samplesheet (single row) -----------
    // Pass the VCF tuple directly — GENERATE_SAMPLESHEET resolves the real path
    GENERATE_SAMPLESHEET(JOINT_GENOTYPE.out.vcf)

    emit:
    samplesheet = GENERATE_SAMPLESHEET.out.csv   // path to CSV
    joint_vcf   = JOINT_GENOTYPE.out.vcf          // tuple(vcf, tbi)
}