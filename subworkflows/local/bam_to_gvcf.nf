include { LINK_LOCAL_BAM          } from '../../modules/local/bam_to_gvcf/link_local_bam'
include { GDC_DOWNLOAD            } from '../../modules/local/bam_to_gvcf/gdc_download'
include { URL_DOWNLOAD            } from '../../modules/local/bam_to_gvcf/url_download'
include { REF_GENOME_RECOGNITION  } from '../../modules/local/bam_to_gvcf/ref_genome_recognition'
include { REALIGN_BWA_MEM2        } from '../../modules/local/bam_to_gvcf/realign_bwa_mem2'
include { COORDINATE_SORT         } from '../../modules/local/bam_to_gvcf/coordinate_sort'
include { PREPARE_COHORT_REF      } from '../../modules/local/bam_to_gvcf/prepare_cohort_reference'
include { FILTER_DBSNP            } from '../../modules/local/bam_to_gvcf/filter_dbsnp'
include { MARK_DUPLICATES         } from '../../modules/local/bam_to_gvcf/mark_duplicates'
include { BASE_RECALIBRATOR       } from '../../modules/local/bam_to_gvcf/base_recalibrator'
include { GATHER_BQSR_REPORTS     } from '../../modules/local/bam_to_gvcf/gather_bqsr_reports'
include { APPLY_BQSR              } from '../../modules/local/bam_to_gvcf/apply_bqsr'
include { GATHER_BAM_FILES        } from '../../modules/local/bam_to_gvcf/gather_bam_files'
include { WGS_METRICS             } from '../../modules/local/bam_to_gvcf/wgs_metrics'
include { FLAGSTAT                } from '../../modules/local/bam_to_gvcf/flagstat'
include { SPLIT_CALLING_REGIONS   } from '../../modules/local/bam_to_gvcf/split_calling_regions'
include { SCATTER_INTERVALS       } from '../../modules/local/bam_to_gvcf/scatter_intervals'
include { HAPLOTYPECALLER         } from '../../modules/local/bam_to_gvcf/haplotypecaller'
include { MERGE_GVCFS             } from '../../modules/local/bam_to_gvcf/merge_gvcfs'
include { JOINT_GENOTYPE          } from '../../modules/local/bam_to_gvcf/joint_genotype'
include { GENERATE_SAMPLESHEET    } from '../../modules/local/bam_to_gvcf/generate_samplesheet'

workflow BAM_TO_GVCF {

    take:
    ch_sample_table   // channel: path to TSV with sampleId (optional) and bamFile columns

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

    // ---- Parse sample table -----------------------------------------------  // *** CHANGED ***
    // Now emits (sampleId, bamFile) tuples.
    // 'sampleId' is read from the TSV 'sampleId' column when present;
    // falls back to the filename portion of bamFile (minus .bam extension).
    ch_all_ids = ch_sample_table
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            if (!row.bamFile) error "Missing 'bamFile' column in row: ${row}"
            def bamFile  = row.bamFile as String
            def sampleId = row.sampleId                                          // null if column absent
                ?: bamFile.tokenize('/').last().replaceAll(/\.bam$/, '')         // filename fallback
            tuple(sampleId, bamFile)
        }

    // ---- Branch: URL vs local BAM path vs remote GDC ID -------------------  // *** CHANGED ***
    //  The closure now destructures the (sampleId, bamFile) tuple.
    //  Check for URLs first so that file(bamFile).exists() is never called on an
    //  HTTP/FTP address (Nextflow's file() *can* resolve remote URIs, which
    //  would short-circuit into the local branch).
    ch_all_ids.branch { sampleId, bamFile ->
        url:    bamFile.startsWith('http://') || bamFile.startsWith('https://') || bamFile.startsWith('ftp://')
        local:  file(bamFile).exists()
        remote: true                                                             // GDC UUIDs
    }.set { sample_branch }

    // -- Local BAMs — resolve BAI companion ---------------------------------  // *** CHANGED ***
    ch_local_input = sample_branch.local.map { sampleId, bam_path ->
        // sampleId now comes from the TSV, not derived from the filename
        def bam = file(bam_path)
        def bai = file("${bam}.bai").exists()
                ? file("${bam}.bai")
                : ( file(bam_path.replaceAll(/\.bam$/, '.bai')).exists()
                    ? file(bam_path.replaceAll(/\.bam$/, '.bai'))
                    : null )
        if (!bai) error "Cannot find BAI for: ${bam_path}"
        tuple(sampleId, bam, bai)
    }

    LINK_LOCAL_BAM(ch_local_input)

    // -- URL BAMs — sampleId comes from the TSV (or filename fallback) ------  // *** CHANGED ***
    // sample_branch.url is already a (sampleId, bam_url) tuple — pass directly.
    // The old ch_url_input mapping that re-derived sampleId from the filename is removed.
    URL_DOWNLOAD(sample_branch.url)

    // -- GDC BAMs -----------------------------------------------------------
    // NOTE: GDC_DOWNLOAD must also be updated to accept a (sampleId, uuid) tuple
    GDC_DOWNLOAD(sample_branch.remote)

    // ---- Unify all three sources ------------------------------------------
    ch_downloaded = LINK_LOCAL_BAM.out.bam
        .mix(URL_DOWNLOAD.out.bam)
        .mix(GDC_DOWNLOAD.out.bam)

    // ---- Reference genome recognition -------------------------------------
    REF_GENOME_RECOGNITION(ch_downloaded, ch_ref_dict)

    REF_GENOME_RECOGNITION.out.result.branch {
        to_realign: it[1] != '0'   // needsRealign == 1 → ('yes'): mismatch found
        keep:       true           // needsRealign == 0 → ('no'): already matches reference
    }.set { branched_bams }

    ch_to_realign = branched_bams.to_realign.map { sid, needsRealign, bam, bai -> tuple(sid, bam, bai) }
    ch_keep       = branched_bams.keep.map       { sid, needsRealign, bam, bai -> tuple(sid, bam, bai) }

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


    // ====================================================================
    //  Dedup + BQSR
    // ====================================================================

    // Step 1: MarkDuplicates — whole genome, per sample (cannot scatter)
    MARK_DUPLICATES(COORDINATE_SORT.out.bam)

    // Step 2: BaseRecalibrator — scattered across 22 autosomes per sample
    ch_bqsr_scatter = MARK_DUPLICATES.out.bam                     // (sid, bam, bai)
        .combine(PREPARE_COHORT_REF.out.reference)                // + (ref, fai, dict)
        .combine(FILTER_DBSNP.out.dbsnp)                          // + (dbsnp, tbi)
        .combine(Channel.of(*CHROMS))                              // × 22 chroms
        .map { sid, bam, bai, ref, ref_fai, ref_dict,
               dbsnp, dbsnp_tbi, chrom ->
            tuple(sid, bam, bai, ref, ref_fai, ref_dict,
                  dbsnp, dbsnp_tbi, chrom)
        }

    BASE_RECALIBRATOR(ch_bqsr_scatter)

    // Step 3: Gather per-chromosome recalibration tables → one table per sample
    BASE_RECALIBRATOR.out.table
        .groupTuple(by: 0)                                        // (sid, [tables])
        .set { ch_grouped_tables }

    GATHER_BQSR_REPORTS(ch_grouped_tables)

    // Step 4: ApplyBQSR — scattered across 22 autosomes per sample
    ch_apply_scatter = MARK_DUPLICATES.out.bam
        .join(GATHER_BQSR_REPORTS.out.table)                      // (sid, bam, bai, recal)
        .combine(PREPARE_COHORT_REF.out.reference)                // + (ref, fai, dict)
        .combine(Channel.of(*CHROMS))                              // × 22 chroms
        .map { sid, bam, bai, recal_table,
               ref, ref_fai, ref_dict, chrom ->
            tuple(sid, bam, bai, ref, ref_fai, ref_dict,
                  recal_table, chrom)
        }

    APPLY_BQSR(ch_apply_scatter)

    // Step 5: Gather per-chromosome BAMs → one final BAM per sample
    APPLY_BQSR.out.bam                                            // (sid, chrom, bam, bai)
        .groupTuple(by: 0)                                        // (sid, [chroms], [bams], [bais])
        .map { sid, chroms, bams, bais -> tuple(sid, bams, bais) }
        .set { ch_gather_bams }

    GATHER_BAM_FILES(ch_gather_bams)

    // ====================================================================


    // ---- Optional QC metrics ----------------------------------------------
    if (params.bam2gvcf_run_metrics != false) {
        ch_metrics_input = GATHER_BAM_FILES.out.bam
            .combine(PREPARE_COHORT_REF.out.reference)
            .map { sid, bam, bai, ref, ref_fai, ref_dict ->
                tuple(sid, bam, bai, ref, ref_fai, ref_dict)
            }
        WGS_METRICS(ch_metrics_input)
        FLAGSTAT(GATHER_BAM_FILES.out.bam)
    }


    // ====================================================================
    //  HaplotypeCaller  (sub-chromosome scatter-gather)
    // ====================================================================

    // ---- Split calling regions by chromosome ------------------------------
    ch_chromosomes = Channel.of(*CHROMS)
    SPLIT_CALLING_REGIONS(ch_chromosomes.combine(ch_calling_regions))

    ch_valid_regions = SPLIT_CALLING_REGIONS.out.bed
        .filter { chr, bed -> bed.size() > 0 }
        .ifEmpty { error "No valid calling regions for any chromosome" }

    // Split each per-chromosome BED into balanced sub-intervals
    SCATTER_INTERVALS(
        ch_valid_regions,                                          // (chrom, bed)
        ch_ref_fasta, ch_ref_fai, ch_ref_dict,
        Channel.value(params.bam2gvcf_hc_scatter_count)
    )

    // Flatten: (chrom, [files]) → (chrom, scatter_idx, file)
    ch_scattered = SCATTER_INTERVALS.out.intervals
        .flatMap { chrom, intervals ->
            def files = intervals instanceof List ? intervals : [intervals]
            files.sort { a, b -> a.name <=> b.name }
                 .withIndex()
                 .collect { f, idx -> tuple(chrom, idx, f) }
        }

    // Build per-sample base channel (now from GATHER_BAM_FILES)
    ch_hc_base = GATHER_BAM_FILES.out.bam
        .combine(PREPARE_COHORT_REF.out.reference)
        .combine(FILTER_DBSNP.out.dbsnp)
        .map { sid, bam, bai, ref, ref_fai, ref_dict,
               dbsnp, dbsnp_tbi ->
            tuple(sid, bam, bai, ref, ref_fai, ref_dict,
                  dbsnp, dbsnp_tbi)
        }

    // Cartesian: sample × (chrom, scatter_idx, interval)
    ch_hc_jobs = ch_hc_base.combine(ch_scattered)
        .map { sid, bam, bai, ref, ref_fai, ref_dict,
               dbsnp, dbsnp_tbi, chrom, scatter_idx, interval ->
            tuple(sid, bam, bai, ref, ref_fai, ref_dict,
                  dbsnp, dbsnp_tbi, chrom, scatter_idx, interval)
        }

    HAPLOTYPECALLER(ch_hc_jobs)

    // Merge sub-interval gVCFs back to per-sample × per-chromosome
    HAPLOTYPECALLER.out.gvcf                                       // (sid, chrom, idx, gvcf, tbi)
        .map { sid, chrom, scatter_idx, gvcf, tbi ->
            tuple(sid, chrom, gvcf, tbi)
        }
        .groupTuple(by: [0, 1])                                   // (sid, chrom, [gvcfs], [tbis])
        .set { ch_merge_gvcfs }

    MERGE_GVCFS(ch_merge_gvcfs)

    MERGE_GVCFS.out.gvcf
        .map { sid, chrom, gvcf, tbi -> tuple(chrom, gvcf, tbi) }
        .groupTuple(by: 0)
        .set { ch_per_chrom }

    // ====================================================================
    //  JOINT_GENOTYPE  (Per-chromosome joint genotyping)
    // ====================================================================

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


    // ====================================================================
    //  GENERATE_SAMPLESHEET
    // ====================================================================

    GENERATE_SAMPLESHEET(ch_all_vcfs)

    emit:
    samplesheet = GENERATE_SAMPLESHEET.out.csv
    joint_vcfs  = JOINT_GENOTYPE.out.vcf
}
