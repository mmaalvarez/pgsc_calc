//
// Check input samplesheet and get read channels
//

include { FORMAT_SCOREFILES } from '../../modules/local/format_scorefiles'


workflow INPUT_CHECK {
    take:
    input_path // file: /path/to/samplesheet.csv  — OR a channel (from BAM/gVCF preprocessing)
    format // csv or JSON
    scorefile // flat list of paths
    chain_files

    main:
    /* all genomic data should be represented as a list of : [[meta], file]

       meta hashmap structure:
        id: experiment label, possibly shared across split genomic files
        is_vcf: boolean, is in variant call format
        is_bfile: boolean, is in PLINK1 fileset format
        is_pfile: boolean, is in PLINK2 fileset format
        chrom: The chromosome associated with the file. If multiple chroms, null.
        n_chrom: Total separate chromosome files per experiment ID
     */

    ch_versions = Channel.empty()
    parsed_input = Channel.empty()

    // ---- Normalise input: string path vs. channel from preprocessing ------
    def is_path_string = (input_path instanceof String || input_path instanceof GString)

    if (is_path_string) {
        input = Channel.fromPath(input_path, checkIfExists: true)
    } else {
        // Already a channel (e.g. from BAM_TO_GVCF / GVCF_TO_JOINT)
        input = input_path
    }

    if (format.equals("csv")) {

        if (is_path_string) {
            // ---- Original eager logic (string path from params.input) ------
            def n_chrom
            n_chrom = file(input_path).countLines() - 1 // ignore header
            parser = new SamplesheetParser(file(input_path), n_chrom, params.target_build)
            input.splitCsv(header:true)
                    .collect()
                    .map { rows -> parser.verifySamplesheet(rows) }
                    .flatten()
                    .map { row -> parser.parseCSVRow(row)}
                    .set { parsed_input }

        } else {
            // ---- Reactive logic (channel from BAM/gVCF preprocessing) ------
            //  Cannot call file() or countLines() eagerly on a channel, so
            //  parse everything inside .map where the concrete Path is available.
            input
                .map { f ->
                    def lines   = f.readLines()
                    def n_chrom = lines.size() - 1          // ignore header
                    def parser  = new SamplesheetParser(f, n_chrom, params.target_build)
                    def header  = lines[0].split(',')*.trim()
                    def rows    = lines.drop(1).collect { line ->
                        def vals = line.split(',', -1)*.trim()
                        [header, vals].transpose().collectEntries()
                    }
                    parser.verifySamplesheet(rows)
                    rows.collect { row -> parser.parseCSVRow(row) }
                }
                .flatMap { it }
                .set { parsed_input }
        }

    } else if (format.equals("json")) {
        // JSON is only used in normal (string-path) mode
        def n_chrom
        n_chrom = file(input_path).countJson()
        parser = new SamplesheetParser(file(input_path), n_chrom, params.target_build)
        input.splitJson()
            .collect()
            .map { jsonarray -> parser.verifySamplesheet(jsonarray)}
            .flatten()
            .map { json -> parser.parseJSON(json)}
            .set { parsed_input }
    }

    parsed_input.branch {
                vcf: it[0].is_vcf
                bfile: it[0].is_bfile
                pfile: it[0].is_pfile
        }
        .set { ch_branched }

    // branch is like a switch statement, so only one bed / bim was being
    // returned
    ch_branched.bfile.multiMap { it ->
        bed: [it[0], it[1][0]]
        bim: [it[0], it[1][1]]
        fam: [it[0], it[1][2]]
    }
        .set { ch_bfiles }

    ch_branched.pfile.multiMap { it ->
        pgen: [it[0], it[1][0]]
        pvar: [it[0], it[1][1]]
        psam: [it[0], it[1][2]]
    }
        .set { ch_pfiles }

    FORMAT_SCOREFILES ( scorefile, chain_files )

    versions = ch_versions.mix(FORMAT_SCOREFILES.out.versions)

    ch_bfiles.bed.mix(ch_pfiles.pgen).dump(tag: 'geno').set { geno }
    ch_bfiles.bim.mix(ch_pfiles.pvar).dump(tag: 'variants').set { variants }
    ch_bfiles.fam.mix(ch_pfiles.psam).dump(tag: 'pheno').set { pheno }
    ch_branched.vcf.dump(tag: 'input').set{vcf}
    FORMAT_SCOREFILES.out.scorefiles.dump(tag: 'input').set{ scorefiles }
    FORMAT_SCOREFILES.out.log_scorefiles.dump(tag: 'input').set{ log_scorefiles }

    emit:
    geno
    variants
    pheno
    vcf
    scorefiles
    log_scorefiles
    versions
}