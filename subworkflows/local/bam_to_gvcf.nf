include { LINK_LOCAL_BAM          } from '../../modules/local/bam_to_gvcf/link_local_bam'
include { GDC_DOWNLOAD            } from '../../modules/local/bam_to_gvcf/gdc_download'
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
include { GATHER_GVCFS            } from '../../modules/local/bam_to_gvcf/gather_gvcfs'
include { GENERATE_SAMPLESHEET    } from '../../modules/local/bam_to_gvcf/generate_samplesheet'

workflow BAM_TO_GVCF {

    take:
    ch_sample_table   // channel: path to TSV with bamFile column

    main:

    // ---- Constants --------------------------------------------------------
    def CHROMS = [
        'chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10',
        'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20',
        'chr21','chr22','chrX','chrY'
    ]
    def chromOrder = CHROMS.withIndex().collectEntries { c, i -> [(c): i] }

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

    // ---- Branch: local BAM path vs remote GDC ID --------------------------
    ch_all_ids.branch {
        local:  file(it).exists()
        remote: true
    }.set { sample_branch }

    // Local BAMs — resolve BAI companion
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
    GDC_DOWNLOAD(sample_branch.remote)

    // ---- Unify both sources -----------------------------------------------
    ch_downloaded = LINK_LOCAL_BAM.out.bam.mix(GDC_DOWNLOAD.out.bam)

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
    ch_first_bam = COORDINATE_SORT.out.bam.first()
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

    // ---- HaplotypeCaller per chromosome -----------------------------------
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

    // ---- Gather per-sample gVCFs across chromosomes -----------------------
    ch_grouped = HAPLOTYPECALLER.out.gvcf
        .groupTuple(by: 0)
        .map { sid, chrom_list, gvcf_list, tbi_list ->
            def zipped = [chrom_list.toList(), gvcf_list.toList(), tbi_list.toList()].transpose()
            zipped.sort { a, b ->
                (chromOrder[a[0].toString()] ?: 9999) <=> (chromOrder[b[0].toString()] ?: 9999)
            }
            tuple(sid, zipped.collect{it[0]}, zipped.collect{it[1]}, zipped.collect{it[2]})
        }

    GATHER_GVCFS(ch_grouped)

    // ---- Generate pgsc_calc-compatible samplesheet ------------------------
    ch_gvcf_names = GATHER_GVCFS.out.gvcf
        .map { sid, gvcf, tbi -> gvcf.name }
        .collect()
        .map { it.sort() }

    GENERATE_SAMPLESHEET(ch_gvcf_names)

    emit:
    samplesheet = GENERATE_SAMPLESHEET.out.csv   // path to CSV
    gvcfs       = GATHER_GVCFS.out.gvcf           // tuple(sampleId, gvcf, tbi)
}