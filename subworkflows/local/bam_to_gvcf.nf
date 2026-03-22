include { LINK_LOCAL_BAM          } from '../../modules/local/bam_to_gvcf/link_local_bam'
include { GDC_DOWNLOAD            } from '../../modules/local/bam_to_gvcf/gdc_download'
include { URL_DOWNLOAD            } from '../../modules/local/bam_to_gvcf/url_download'          // *** NEW ***
include { REF_GENOME_RECOGNITION  } from '../../modules/local/bam_to_gvcf/ref_genome_recognition'
include { REALIGN_BWA_MEM2        } from '../../modules/local/bam_to_gvcf/realign_bwa_mem2'
include { COORDINATE_SORT         } from '../../modules/local/bam_to_gvcf/coordinate_sort'
include { PREPARE_COHORT_REF      } from '../../modules/local/bam_to_gvcf/prepare_cohort_reference'
include { FILTER_DBSNP            } from '../../modules/local/bam_to_gvcf/filter_dbsnp'
include { DEDUP_BQSR              } from '../../modules/local/bam_to_gvcf/dedup_bqsr'
include { WGS_METRICS             } from '../../modules/local/bam_to_gvcf/wgs_metrics'
include { FLAGSTAT                } from '../../modules/local/bam_to_gvcf/flagstat'
include { SPLIT_CALLING_REGIONS   } from '../../modules/local/bam_to_gvcf/split_calling_regions'
include { HAPLOTYPECALLER         } from '../../modules/local/bam_to_gvcf/haplotypecaller'
include { JOINT_GENOTYPE          } from '../../modules/local/bam_to_gvcf/joint_genotype'
include { GENERATE_SAMPLESHEET    } from '../../modules/local/bam_to_gvcf/generate_samplesheet'

workflow BAM_TO_GVCF {

    take:
    ch_sample_table   // channel: path to TSV with bamFile column

    main:

    // ---- Constants --------------------------------------------------------
    // ONLY AUTOSOMES
    def CHROMS = [
        'chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10',
        'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20',
        'chr21','chr22'
    ]

    // ---- Reference files from params --------------------------------------
    ch_ref_fasta   = Channel.value(file(params.bam2gvcf_fasta, checkIfExists: true))
    ch_ref_fai     = Channel.value(file("${params.bam2gvcf_fasta}.fai", checkIfExists: true))
    ch_ref_dict    = Channel.value(file(params.bam2gvcf_dict, checkIfExists: true))

    ch_bwa_pac     = Channel.value(file("${params.bam2gvcf_fasta}.pac", checkIfExists: true))
    ch_bwa_ann     = Channel.value(file("${params.bam2gvcf_fasta}.ann", checkIfExists: true))
    ch_bwa_amb     = Channel.value(file("${params.bam2gvcf_fasta}.amb", checkIfExists: true))
    ch_bwa_0123    = Channel.value(file("${params.bam2gvcf_fasta}.0123", checkIfExists: true))
    ch_bwa_bwt2bit = Channel.value(file("${params.bam2gvcf_fasta}.bwt.2bit.64", checkIfExists: true))

    ch_dbsnp       = Channel.value(file(params.bam2gvcf_dbsnp, checkIfExists: true))
    ch_dbsnp_tbi   = Channel.value(file("${params.bam2gvcf_dbsnp}.tbi", checkIfExists: true))

    ch_calling_regions = Channel.value(file(params.bam2gvcf_calling_regions, checkIfExists: true))
    ch_contig_map      = Channel.value(file(params.bam2gvcf_contig_map, checkIfExists: true))

    // ---- Parse sample table -----------------------------------------------
    ch_all_ids = ch_sample_table
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            if (!row.bamFile) error "Missing 'bamFile' column in row: ${row}"
            row.bamFile as String
        }

    // ---- Branch: URL vs local BAM path vs remote GDC ID -------------------  // *** CHANGED ***
    //  Check for URLs first so that file(it).exists() is never called on an
    //  HTTP/FTP address (Nextflow's file() *can* resolve remote URIs, which
    //  would short-circuit into the local branch).
    ch_all_ids.branch {
        url:    it.startsWith('http://') || it.startsWith('https://') || it.startsWith('ftp://')   // *** NEW ***
        local:  file(it).exists()
        remote: true                                                                                // GDC UUIDs
    }.set { sample_branch }

    // -- Local BAMs — resolve BAI companion ---------------------------------
    ch_local_input = sample_branch.local.map { bam_path ->
        def bam = file(bam_path)
        def sampleId = bam.baseName
        def bai = file("${bam}.bai").exists()
                ? file("${bam}.bai")
                : ( file(bam_path.replaceAll(/\.bam$/, '.bai')).exists()
                    ? file(bam_path.replaceAll(/\.bam$/, '.bai'))
                    : null )
        if (!bai) error "Cannot find BAI for: ${bam_path}"
        tuple(sampleId, bam, bai)
    }

    LINK_LOCAL_BAM(ch_local_input)

    // -- URL BAMs — derive sampleId from the filename in the URL ------------  // *** NEW ***
    ch_url_input = sample_branch.url.map { bam_url ->
        def filename = bam_url.tokenize('/').last()                  // e.g. HG00096.mapped.ILLUMINA.bwa.GBR.low_coverage.20120522.bam
        def sampleId = filename.replaceAll(/\.bam$/, '')
        tuple(sampleId, bam_url)
    }

    URL_DOWNLOAD(ch_url_input)                                                                      // *** NEW ***

    // -- GDC BAMs -----------------------------------------------------------
    GDC_DOWNLOAD(sample_branch.remote)

    // ---- Unify all three sources ------------------------------------------  // *** CHANGED ***
    ch_downloaded = LINK_LOCAL_BAM.out.bam
        .mix(URL_DOWNLOAD.out.bam)                                                                  // *** NEW ***
        .mix(GDC_DOWNLOAD.out.bam)

    // ---- Reference genome recognition -------------------------------------
    REF_GENOME_RECOGNITION(ch_downloaded, ch_ref_dict)

    REF_GENOME_RECOGNITION.out.result.branch {
        to_realign: it[1] != '0'
        keep:       true
    }.set { branched_bams }

    ch_to_realign = branched_bams.to_realign.map { sid, flag, bam, bai -> tuple(sid, bam, bai) }
    ch_keep       = branched_bams.keep.map       { sid, flag, bam, bai -> tuple(sid, bam, bai) }

    // ---- Realignment ------------------------------------------------------
    REALIGN_BWA_MEM2(
        ch_to_realign,
        ch_ref_fasta,
        ch_bwa_pac, ch_bwa_ann, ch_bwa_amb, ch_bwa_0123, ch_bwa_bwt2bit
    )

    ch_post_align = REALIGN_BWA_MEM2.out.bam.mix(ch_keep)

    // ---- Coordinate sort --------------------------------------------------
    COORDINATE_SORT(ch_post_align)

    // ---- Cohort reference from first BAM ----------------------------------
    ch_first_bam = COORDINATE_SORT.out.bam
        .toSortedList { a, b -> a[0] <=> b[0] }  // sort by sampleId
        .map { it.first() }                      // always picks alphabetically first sample

    PREPARE_COHORT_REF(ch_first_bam, ch_ref_fasta, ch_contig_map)

    // ---- Filter dbSNP for cohort ------------------------------------------
    FILTER_DBSNP(PREPARE_COHORT_REF.out.reference, ch_dbsnp, ch_dbsnp_tbi)

    // ---- Dedup + BQSR -----------------------------------------------------
    ch_dedup_input = COORDINATE_SORT.out.bam
        .combine(PREPARE_COHORT_REF.out.reference)
        .combine(FILTER_DBSNP.out.dbsnp)
        .map { sid, bam, bai, ref, ref_fai, ref_dict, dbsnp, dbsnp_tbi ->
            tuple(sid, bam, bai, ref, ref_fai, ref_dict, dbsnp, dbsnp_tbi)
        }

    DEDUP_BQSR(ch_dedup_input)

    // ---- Optional QC metrics ----------------------------------------------
    if (params.bam2gvcf_run_metrics != false) {
        ch_metrics_input = DEDUP_BQSR.out.bam
            .combine(PREPARE_COHORT_REF.out.reference)
            .map { sid, bam, bai, ref, ref_fai, ref_dict ->
                tuple(sid, bam, bai, ref, ref_fai, ref_dict)
            }
        WGS_METRICS(ch_metrics_input)
        FLAGSTAT(DEDUP_BQSR.out.bam)
    }

    // ---- Split calling regions by chromosome ------------------------------
    ch_chromosomes = Channel.of(*CHROMS)
    SPLIT_CALLING_REGIONS(ch_chromosomes.combine(ch_calling_regions))

    ch_valid_regions = SPLIT_CALLING_REGIONS.out.bed
        .filter { chr, bed -> bed.size() > 0 }
        .ifEmpty { error "No valid calling regions for any chromosome" }

    // ---- HaplotypeCaller per sample × chromosome --------------------------
    ch_hc_base = DEDUP_BQSR.out.bam
        .combine(PREPARE_COHORT_REF.out.reference)
        .combine(FILTER_DBSNP.out.dbsnp)
        .map { sid, bam, bai, ref, ref_fai, ref_dict, dbsnp, dbsnp_tbi ->
            tuple(sid, bam, bai, ref, ref_fai, ref_dict, dbsnp, dbsnp_tbi)
        }

    ch_hc_jobs = ch_hc_base.combine(ch_valid_regions)
        .map { sid, bam, bai, ref, ref_fai, ref_dict, dbsnp, dbsnp_tbi, chrom, bed ->
            tuple(sid, bam, bai, ref, ref_fai, ref_dict, dbsnp, dbsnp_tbi, chrom, bed)
        }

    HAPLOTYPECALLER(ch_hc_jobs)

    // ---- Group HaplotypeCaller output BY CHROMOSOME (across all samples) --
    HAPLOTYPECALLER.out.gvcf
        .map { sid, chrom, gvcf, tbi -> tuple(chrom, gvcf, tbi) }
        .groupTuple(by: 0)
        .set { ch_per_chrom }

    // ---- Per-chromosome joint genotyping ----------------------------------
    ch_per_chrom
        .combine(PREPARE_COHORT_REF.out.reference)
        .combine(FILTER_DBSNP.out.dbsnp)
        .map { chrom, gvcfs, tbis, ref, ref_fai, ref_dict, dbsnp, dbsnp_tbi ->
            tuple(chrom, gvcfs, tbis, ref, ref_fai, ref_dict, dbsnp, dbsnp_tbi)
        }
        .set { ch_jg_input }

    JOINT_GENOTYPE(ch_jg_input)

    // ---- Generate pgsc_calc samplesheet (one row per chromosome) ----------
    JOINT_GENOTYPE.out.vcf
        .map { chrom, vcf, tbi -> vcf }
        .collect()
        .set { ch_all_vcfs }

    GENERATE_SAMPLESHEET(ch_all_vcfs)

    emit:
    samplesheet = GENERATE_SAMPLESHEET.out.csv
    joint_vcfs  = JOINT_GENOTYPE.out.vcf
}
