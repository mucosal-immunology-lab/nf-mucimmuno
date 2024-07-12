// IMPORT PROCESSES

include { FASTQC } from './modules/fastqc/main.nf'
include { TRIMGALORE } from './modules/trimgalore/main.nf'
include { MULTIQC as MULTIQC_PRETRIM } from './modules/multiqc/main.nf'
include { MULTIQC as MULTIQC_POSTTRIM } from './modules/multiqc/main.nf'
include { STARSOLO } from './modules/starsolo/main.nf'
include { COLLECT_EXPORT_FILES } from './modules/collect_export_files/main.nf'


workflow {
    // Create input channel from the samplesheet provided through params.samples_csv
    Channel.fromPath(params.samples_csv, checkIfExists:true)
    | splitCsv(header:true)
    | map { row ->
        meta = row.subMap('id') + [CLS: row.CLS] + [GenomeIndex: row.GenomeIndex]
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

    // =====================================================
    // R U N    S T A R s o l o    O N    R A W    R E A D S
    // =====================================================
    starsolo_counts = Channel.empty()
    starsolo_log_final = Channel.empty()
    starsolo_summary = Channel.empty()
    ch_starsolo = STARSOLO(trimgalore_reads)
    starsolo_counts = ch_starsolo.counts
    starsolo_log_final = ch_starsolo.log_final
    starsolo_summary = ch_starsolo.summary

    // ==================================================
    // C O L L E C T    F I L E S    F O R    E X P O R T
    // ==================================================

    // Collect the STARsolo output folders into a single list
    starsolo_counts
    | map { it[1] }
    | collect
    | set { export_files }

    // Collect the files for export
    COLLECT_EXPORT_FILES(export_files)
}