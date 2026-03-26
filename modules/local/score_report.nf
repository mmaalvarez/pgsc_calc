process SCORE_REPORT {
    // first element of tag must be sampleset
    tag "$meta.id" 

    label 'process_high_memory'
    label 'error_retry'
    label 'report'

    conda "${task.ext.conda}"

    container "${ workflow.containerEngine == 'singularity' &&
        !task.ext.singularity_pull_docker_container ?
        "${task.ext.singularity}${task.ext.singularity_version}" :
        "${task.ext.docker}${task.ext.docker_version}" }"

    beforeScript = {
        def clean = task.attempt > 1 ?
            "CONDA_PKGS_DIRS=\$(cd ../.. && pwd)/conda/pkgs conda clean --all -y && " : ""
        clean + """
        # Derive the active conda env root from the activated Rscript binary
        _CP=\$(dirname \$(dirname \$(which Rscript)))
        # Positively set R_HOME and all library paths to the conda env — no unsetting, no guessing
        export R_HOME=\$_CP/lib/R
        export R_LIBS=\$_CP/lib/R/library
        export R_LIBS_SITE=\$_CP/lib/R/library
        unset R_LIBS_USER
        # Suppress ~/.Renviron and ~/.Rprofile so they can't re-inject a bad R_LIBS_USER
        export R_ENVIRON_USER=/dev/null
        export R_PROFILE_USER=/dev/null
        # Tell Quarto explicitly which R binary to use
        export QUARTO_R=\$(which R)
        # Debug: will appear in .command.err so you can confirm the paths
        echo "DEBUG R_HOME=\$R_HOME QUARTO_R=\$QUARTO_R" >&2
        Rscript --no-environ --no-site-file --no-init-file -e ".libPaths()" >&2
        """
    }
    
    input:
    tuple val(meta), path(scorefile), path(score_log), path(match_summary), path(ancestry)
    path intersect_count
    val reference_panel_name
    path(report_path, arity: '4') // 4 files expected: report, css, background image x2

    output:
    // includeInputs to correctly use $meta.id in publishDir path
    // ancestry results are optional also
    path "*.txt.gz", includeInputs: true
    path "*.json.gz", includeInputs: true, optional: true
    // for testing ancestry workflow
    path "pop_summary.csv", optional: true
    // normal outputs
    path "report.html", emit: report
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    run_ancestry = params.run_ancestry ? true : false
    """
    export TMPDIR=\$PWD # tmpdir must always be writable for quarto
    echo $workflow.commandLine > command.txt
    
    echo "keep_multiallelic: $params.keep_multiallelic" > params.txt
    echo "keep_ambiguous   : $params.keep_ambiguous"    >> params.txt
    echo "min_overlap      : $params.min_overlap"       >> params.txt
    
    quarto render report.qmd -M "self-contained:true" \
        -P score_path:$scorefile \
        -P sampleset:$meta.id \
        -P run_ancestry:$run_ancestry \
        -P reference_panel_name:$reference_panel_name \
        -P version:$workflow.manifest.version \
        -o report.html

    cat <<-END_VERSIONS > versions.yml
    ${task.process.tokenize(':').last()}:
        R: \$(echo \$(R --version 2>&1) | head -n 1 | cut -f 3 -d ' ')
    END_VERSIONS
    """
}

