process REALIGN_BWA_MEM2 {
    tag "${sampleId}"

    conda "bioconda::bwa-mem2=2.2.1 bioconda::samtools=1.17 bioconda::pysam=0.22.0"

    def cachedir = params.genotypes_cache ? file(params.genotypes_cache) : workDir
    storeDir cachedir / "bam_to_gvcf" / "bam_realigned"

    input:
    tuple val(sampleId), path(bam_file), path(bai_file)
    path(reference)
    path(reference_pac)
    path(reference_ann)
    path(reference_amb)
    path(reference_0123)
    path(reference_bwt2bit)

    output:
    tuple val(sampleId),
          path("${sampleId}.realigned.sorted.bam"),
          path("${sampleId}.realigned.sorted.bam.bai"), emit: bam

    shell:
    '''
    set -euo pipefail

    # ── Pre-validate input BAM ──────────────────────────────────────────
    echo "[INFO] Checking input BAM integrity..."
    samtools quickcheck !{bam_file} || {
        echo "ERROR: input BAM !{bam_file} fails integrity check (truncated/corrupt?)" >&2
        exit 1
    }

    # ── Setup ───────────────────────────────────────────────────────────
    tmp_dir_collate=$(mktemp -d --tmpdir=. collate.XXXX)
    tmp_dir_sort=$(mktemp -d --tmpdir=. sort.XXXX)

    bwa_threads=$(( !{task.cpus} - 9 ))
    if [ $bwa_threads -lt 1 ]; then bwa_threads=1; fi

    readGroupInfo=$(retrieve_read_group_info.py !{bam_file})

    # ── Named pipe so we can capture the extraction exit code ───────────
    mkfifo interleaved.fifo

    # BAM → interleaved FASTQ runs in a background subshell.
    # pipefail inside the subshell ensures that if samtools collate fails,
    # the whole subshell fails (not just samtools fastq).
    (
        set -euo pipefail
        samtools collate -f -O !{bam_file} -@ 2 "$tmp_dir_collate" | \
            samtools fastq -@ 2 - > interleaved.fifo
    ) &
    pid_extract=$!

    # ── Alignment pipeline reads from the FIFO ─────────────────────────
    bwa-mem2 mem -t "$bwa_threads" !{reference} -p interleaved.fifo | \
        samtools addreplacerg -r "@RG\\tID:!{sampleId}\\tSM:!{sampleId}\\t$readGroupInfo" - | \
        samtools view -Sb - | \
        samtools sort -o !{sampleId}.realigned.sorted.bam -@ 5 -T "$tmp_dir_sort" -

    # ── Check whether the extraction actually succeeded ─────────────────
    wait $pid_extract || {
        echo "ERROR: BAM-to-FASTQ extraction failed for !{bam_file} (corrupt/truncated input?)" >&2
        rm -f !{sampleId}.realigned.sorted.bam
        rm -rf "$tmp_dir_collate" "$tmp_dir_sort" interleaved.fifo
        exit 1
    }

    # ── Index and clean up ──────────────────────────────────────────────
    samtools index !{sampleId}.realigned.sorted.bam

    rm -rf "$tmp_dir_collate" "$tmp_dir_sort" interleaved.fifo
    '''
}