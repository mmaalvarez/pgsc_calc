process MARK_DUPLICATES {
    tag "${sampleId}"

    conda "bioconda::gatk4=4.5.0.0"

    input:
    tuple val(sampleId), path(coord_sorted_bam), path(coord_sorted_bai)

    output:
    tuple val(sampleId),
          path("${sampleId}.dedup.bam"),
          path("${sampleId}.dedup.bam.bai"), emit: bam
    path("${sampleId}.marked_dup_metrics.txt"),  emit: metrics

    shell:
    '''
    set -euo pipefail
    tmp_dir=$(mktemp -d --tmpdir=. md.XXXX)

    gatk --java-options "-Xmx!{Math.max((task.memory.toGiga() as int)-2,1)}G" MarkDuplicates \
      -I !{coord_sorted_bam} \
      -O !{sampleId}.dedup.bam \
      -M !{sampleId}.marked_dup_metrics.txt \
      --REMOVE_DUPLICATES true \
      --TMP_DIR $tmp_dir

    gatk BuildBamIndex \
      -I !{sampleId}.dedup.bam \
      -O !{sampleId}.dedup.bam.bai

    rm -rf $tmp_dir
    '''
}