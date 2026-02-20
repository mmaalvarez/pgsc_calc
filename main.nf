#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    mmaalvarez/pgsc_calc
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/mmaalvarez/pgsc_calc
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsHelp } from 'plugin/nf-schema'

// Print help message if needed
if (params.help) {
    log.info paramsHelp("nextflow run mmaalvarez/pgsc_calc --input input_file.csv")
    log.info "See https://pgsc-calc.readthedocs.io/en/latest/getting-started.html for more help"
    exit 0
}

WorkflowMain.initialise(workflow, params, log, args)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PGSCCALC }      from './workflows/pgsc_calc'
include { BAM_TO_GVCF }   from './subworkflows/local/bam_to_gvcf'
include { GVCF_TO_JOINT }  from './subworkflows/local/gvcf_to_joint'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HELPER: detect input type from the header line
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def isBamSampleTable(input_path) {
    def first_line = file(input_path).newReader().readLine()
    return first_line.contains('bamFile')
}

def isGvcfSampleTable(input_path) {
    def first_line = file(input_path).newReader().readLine()
    return first_line.contains('gvcfFile')
}

//
// WORKFLOW: Run main pgscatalog/pgsccalc analysis pipeline
//
workflow PGSCATALOG_PGSCCALC {

    if (isBamSampleTable(params.input)) {
        // ---- Mode 1: BAM table → realign/call → joint genotype → pgsc_calc ----
        log.info "Detected BAM sample table — running BAM-to-gVCF preprocessing"

        BAM_TO_GVCF(Channel.fromPath(params.input))
        PGSCCALC(BAM_TO_GVCF.out.samplesheet)

    } else if (isGvcfSampleTable(params.input)) {
        // ---- Mode 2: gVCF table → joint genotype → pgsc_calc ----
        log.info "Detected gVCF sample table — running joint genotyping only"

        GVCF_TO_JOINT(Channel.fromPath(params.input))
        PGSCCALC(GVCF_TO_JOINT.out.samplesheet)

    } else {
        // ---- Mode 3: Normal pgsc_calc samplesheet (VCF / pfile / bfile) ----
        PGSCCALC(params.input)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    PGSCATALOG_PGSCCALC ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
    |\__/,|   (`\
  _.|o o  |_   ) )
-(((---(((--------
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/