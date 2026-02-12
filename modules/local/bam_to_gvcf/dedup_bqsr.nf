process DEDUP_BQSR {
    tag "${sampleId}"

    input:
    tuple val(sampleId), path(coord_sorted_bam), path(coord_sorted_bai),
          path(reference), path(reference_fai), path(reference_dict),
          path(dbsnp), path(dbsnp_index)

    output:
    tuple val(sampleId),
          path("${sampleId}_final.bam"),
          path("${sampleId}_final.bam.bai"), emit: bam

    shell:
    '''
    set -euo pipefail
    tmp_dir_md=$(mktemp -d --tmpdir=. md.XXXX)

    gatk --java-options "-Xmx!{Math.max((task.memory.toGiga() as int)-2,1)}G" MarkDuplicates \
      -I !{coord_sorted_bam} \
      -O !{sampleId}.dedup.bam \
      -M !{sampleId}.marked_dup_metrics.txt \
      --REMOVE_DUPLICATES true \
      --TMP_DIR $tmp_dir_md

    gatk --java-options "-Xmx!{Math.max((task.memory.toGiga() as int)-2,1)}G" BaseRecalibrator \
      -R !{reference} \
      -I !{sampleId}.dedup.bam \
      --known-sites !{dbsnp} \
      -O recal_data.table

    gatk --java-options "-Xmx!{Math.max((task.memory.toGiga() as int)-2,1)}G" ApplyBQSR \
      -R !{reference} \
      -I !{sampleId}.dedup.bam \
      --bqsr-recal-file recal_data.table \
      -O !{sampleId}_final.bam

    gatk BuildBamIndex -I !{sampleId}_final.bam -O !{sampleId}_final.bam.bai
    rm -rf $tmp_dir_md
    '''
}