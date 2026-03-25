process LINK_LOCAL_BAM {
    tag "${sampleId}"

    input:
    tuple val(sampleId), path(bam_file), path(bai_file)

    output:
    // Reference the input path variables directly — never use globs here.
    // This guarantees each emission carries exactly one BAM/BAI regardless
    // of how many files Nextflow stages into the work directory.
    tuple val(sampleId), path(bam_file), path(bai_file), emit: bam

    script:
    """
    set -euo pipefail
    # Staging/symlinking is handled entirely by Nextflow.
    # Nothing else is needed; ls is kept only for debugging.
    ls -la "${bam_file}" "${bai_file}"
    """
}