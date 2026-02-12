process REALIGN_BWA_MEM2 {
    tag "${sampleId}"

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

    tmp_dir_collate=$(mktemp -d --tmpdir=. collate.XXXX)
    tmp_dir_sort=$(mktemp -d --tmpdir=. sort.XXXX)

    bwa_threads=$(( !{task.cpus} - 9 ))
    if [ $bwa_threads -lt 1 ]; then bwa_threads=1; fi

    readGroupInfo=$(retrieve_read_group_info.py !{bam_file})

    bwa-mem2 mem -t $bwa_threads !{reference} -p <( \
        samtools fastq -@ 2 <(samtools collate -f -O !{bam_file} -@ 2 $tmp_dir_collate) \
      ) | \
      samtools addreplacerg -r "@RG\\tID:!{sampleId}\\tSM:!{sampleId}\\t$readGroupInfo" - | \
      samtools view -Sb - | \
      samtools sort -o !{sampleId}.realigned.sorted.bam -@ 5 -T $tmp_dir_sort -

    samtools index !{sampleId}.realigned.sorted.bam

    rm -rf $tmp_dir_collate $tmp_dir_sort
    '''
}