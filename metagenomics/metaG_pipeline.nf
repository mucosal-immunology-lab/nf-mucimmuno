// IMPORT PROCESSES

include { FASTQC } from './modules/fastqc/main.nf'
include { TRIMGALORE } from './modules/trimgalore/main.nf'
include { MULTIQC as MULTIQC_PRETRIM } from './modules/multiqc/main.nf'
include { MULTIQC as MULTIQC_POSTTRIM } from './modules/multiqc/main.nf'
include { KOMPLEXITY_FILTER } from './modules/komplexity/main.nf'
include { PREPARE_HOST_GENOME; HOST_DECONTAMINATE } from './modules/host_decontaminate/main.nf'
include { PREPARE_KRAKEN2_DB; CLASSIFY_KRAKEN2; MERGE_KRAKEN2_REPORTS } from './modules/kraken2/main.nf'
include { PREPARE_BRACKEN_DB; RUN_BRACKEN_CORRECTION; MERGE_BRACKEN_REPORTS } from './modules/bracken/main.nf'

workflow {
    // Create input channel from the samplesheet provided through params.samples_csv
    Channel.fromPath(params.samples_csv, checkIfExists:true)
    | splitCsv(header:true)
    | map { row ->
        meta = row.subMap('id')
        single_end = !row.fastq_2
        files = single_end ? [file(row.fastq_1, checkIfExists:true)] : [file(row.fastq_1, checkIfExists:true), file(row.fastq_2, checkIfExists:true)]
        [meta + [single_end: single_end], files]
    }
    | set { ch_fastq }

    // ===========================================================
    // R U N    F A S T Q C    O N    P R E - T R I M    R E A D S
    // ===========================================================

    // Run FASTQC module from NF-CORE
    fastqc_html_pretrim = Channel.empty()
    fastqc_zip_pretrim = Channel.empty()
    ch_fastqc_pretrim = FASTQC(ch_fastq)
    fastqc_html_pretrim = ch_fastqc_pretrim.html
    fastqc_zip_pretrim = ch_fastqc_pretrim.zip

    // Collect pretrim HTML files into a single list
    fastqc_html_pretrim
    | map { it[1] }
    | collect
    | set { html_pretrim }

    fastqc_zip_pretrim
    | map { it[1] }
    | collect
    | set { zip_pretrim }

    html_pretrim
    | combine(zip_pretrim)
    | collect
    | set { pretrim_reports }

    // Run MultiQC on the pre-trim report
    multiqc_report_pretrim = MULTIQC_PRETRIM(pretrim_reports, 'pretrim')

    // =========================================================
    // R U N    T R I M G A L O R E    O N    R A W    R E A D S
    // =========================================================

    // Run TrimGalore on the raw reads
    trimgalore_reads = Channel.empty()
    trimgalore_html = Channel.empty()
    trimgalore_zip = Channel.empty()
    trimgalore_log = Channel.empty()
    ch_trimgalore = TRIMGALORE(ch_fastq)
    trimgalore_reads = ch_trimgalore.reads
    trimgalore_html = ch_trimgalore.html
    trimgalore_zip = ch_trimgalore.zip
    trimgalore_log = ch_trimgalore.log

    // Collect post-trim HTML files into a single list
    trimgalore_html
    | map { it[1] }
    | collect
    | set { html_posttrim }

    trimgalore_zip
    | map { it[1] }
    | collect
    | set { zip_posttrim }

    html_posttrim
    | combine(zip_posttrim)
    | collect
    | set { posttrim_reports }

    // Run MultiQC on the pre-trim report
    multiqc_report_posttrim = MULTIQC_POSTTRIM(posttrim_reports, 'posttrim')

    // =================================================================
    // R U N    K O M P L E X I T Y    O N    T R I M M E D    R E A D S
    // =================================================================

    // Run Komplexity on the trimmed reads
    komplexity_reads = Channel.empty()
    ch_komplexity = KOMPLEXITY_FILTER(trimgalore_reads)
    komplexity_reads = ch_komplexity.reads

    // ================================================
    // R U N    H O S T   D E C O N T A M I N A T I O N
    // ================================================

    /*
      Determine the host index to use:
      - If `params.decontaminate.hostIndex` is defined in nextflow.config and the file exists,
        then use that as the host index.
      - Otherwise, run the PREPARE_HOST_GENOME process to create the host genome
        and build the Bowtie2 index.
    */
    def hostIndex_ch
    if (params.decontaminate.hostIndex && file(params.decontaminate.hostIndex + ".1.bt2").exists()) {
        log.info "Using provided host index from config: ${params.decontaminate.hostIndex}"
        hostIndex_ch = Channel.value(params.decontaminate.hostIndex)
    } else {
        log.info "No valid host index provided. Running PREPARE_HOST_GENOME process."
        // Launch the process (its output channel emits the index prefix files)
        host_genome_ch = PREPARE_HOST_GENOME()

        // Collect the output from the process and dynamically determine the prefix
        hostIndex_ch = host_genome_ch.bowtie2_index.map { indexFiles ->
            // Extract the common prefix from the Bowtie2 index files
            def prefix = indexFiles[0].toString().replaceAll(/\.1\.bt2$/, '')
            log.info "Generated host index prefix: ${prefix}"
            return prefix
        }
    }

    // Run host decontamination using bowtie2
    host_decontam_reads = Channel.empty()
    ch_host_decontam = HOST_DECONTAMINATE(komplexity_reads, hostIndex_ch)
    host_decontam_reads = ch_host_decontam.reads

    // ================================================================
    // R U N    K R A K E N 2   T A X O N O M I C   A S S I G N M E N T
    // ================================================================

    /*
      Determine the Kraken2 database to use:
      - If `params.taxonomy.kraken2_db` is defined in nextflow.config and the file exists,
        then use that as the kraken2_db.
      - Otherwise, run the PREPARE_KRAKEN2_DB process to create the appropriate
        Kraken2 database.
    */
    def kraken2_db_ch
    if (params.taxonomy.kraken2_db && file(params.taxonomy.kraken2_db).isDirectory()) {
        def kraken2_db_dir = file(params.taxonomy.kraken2_db)
        def requiredFiles = ['hash.k2d', 'database.kraken']
        def missingFiles = requiredFiles.findAll { fileName ->
            !kraken2_db_dir.resolve(fileName).exists()
        }

        if (missingFiles.isEmpty()) {
            log.info "Using provided Kraken2 database from config: ${params.taxonomy.kraken2_db}"
            kraken2_db_ch = Channel.value(params.taxonomy.kraken2_db)
        } else {
            throw new IllegalArgumentException("Provided Kraken2 database directory is missing required files: ${missingFiles.join(', ')}")
        }
    } else {
        log.info "No valid Kraken2 database provided. Running PREPARE_KRAKEN2_DB process."
        // Launch the process (its output channel emits the database path)
        kraken2_db_path = PREPARE_KRAKEN2_DB()

        // Collect the output from the process
        kraken2_db_ch = kraken2_db_path.kraken2_db
    }

    // Run Kraken2 classification
    kraken2_classification = Channel.empty()
    kraken2_report = Channel.empty()
    ch_kraken2 = CLASSIFY_KRAKEN2(host_decontam_reads, kraken2_db_ch)
    kraken2_classification = ch_kraken2.kraken
    kraken2_report = ch_kraken2.report

    // Extract just the report file paths, and collect them
    kraken2_report_paths = kraken2_report.map { meta, rpt -> rpt }
    kraken2_report_paths
        .collect()
        .set { kraken_reports }
        
    // Extract only the IDs, then collect them into a single list
    kraken2_sample_ids = kraken2_report
        .map { meta, rpt -> meta.id }
        .collect()

    // Merge kraken reports together
    kraken2_merged_report = Channel.empty()
    ch_kraken2_merge = MERGE_KRAKEN2_REPORTS(kraken_reports, kraken2_sample_ids)
    kraken2_merged_report = ch_kraken2_merge.combined_report

    // ================================================================
    // R U N    B R A C K E N   T A X O N O M I C   A S S I G N M E N T
    // ================================================================

    def bracken_db_ch
    if ( params.taxonomy.kraken2_db && file(params.taxonomy.kraken2_db).isDirectory() ) {
        def k2db = file(params.taxonomy.kraken2_db)
        // Bracken outputs that must reside in the Kraken2 DB folder
        def requiredFiles = [
            "database${params.taxonomy.kmer_length}mers.kmer_distrib",
            "database${params.taxonomy.kmer_length}mers.kraken"
        ]
        def missing = requiredFiles.findAll { fname -> !k2db.resolve(fname).exists() }

        if ( missing.isEmpty() ) {
            log.info "Using existing Bracken outputs in Kraken2 DB: ${params.taxonomy.kraken2_db}"
            // Point bracken_db_ch to the same directory, downstream processes will load from here
            bracken_db_ch = Channel.value(params.taxonomy.kraken2_db)
        } else {
            log.info "Bracken outputs not found in Kraken2 DB (missing: ${missing.join(', ')}). Running PREPARE_BRACKEN_DB."
            def prep = PREPARE_BRACKEN_DB( kraken2_db: params.taxonomy.kraken2_db )
            bracken_db_ch = prep.bracken_db
        }
    } else {
        log.info "No Kraken2 DB provided or directory not found. Building Bracken DB after Kraken2."
        // Ensure we have kraken2_db_ch defined earlier for the Kraken2 DB
        def prep = PREPARE_BRACKEN_DB( kraken2_db: kraken2_db_ch )
        bracken_db_ch = prep.bracken_db
    }

    // Run Bracken correction
    bracken_report = Channel.empty()
    ch_bracken = RUN_BRACKEN_CORRECTION(kraken2_report)
    bracken_report = ch_bracken.bracken_report

    // Extract just the report file paths, and collect them
    bracken_report_paths = bracken_report.map { meta, rpt -> rpt }
    bracken_report_paths
        .collect()
        .set { bracken_reports }
        
    // Extract only the IDs, then collect them into a single list
    bracken_sample_ids = bracken_report
        .map { meta, rpt -> meta.id }
        .collect()

    // Merge bracken reports together
    bracken_merged_report = Channel.empty()
    ch_bracken_merge = MERGE_BRACKEN_REPORTS(bracken_reports, bracken_sample_ids)
    bracken_merged_report = ch_bracken_merge.combined_report
}