process COORDINATE_SORT {
    tag "${sampleId}"

    input:
    tuple val(sampleId), path(in_bam), path(in_bai)

    output:
    tuple val(sampleId),
          path("${sampleId}.coord.sorted.bam"),
          path("${sampleId}.coord.sorted.bam.bai"), emit: bam

    shell:
    '''
    set -euo pipefail

    so=$(samtools view -H !{in_bam} \
         | awk -F'\\t' '$1=="@HD"{for(i=1;i<=NF;i++) if($i ~ /^SO:/){sub(/^SO:/,"",$i); print $i}}' \
         | head -n 1 || true)

    if [ "$so" = "coordinate" ]; then
      ln -s !{in_bam} !{sampleId}.coord.sorted.bam
      ln -s !{in_bai} !{sampleId}.coord.sorted.bam.bai
    else
      samtools sort -@ !{task.cpus} -o !{sampleId}.coord.sorted.bam !{in_bam}
      samtools index -@ !{task.cpus} !{sampleId}.coord.sorted.bam
    fi
    '''
}