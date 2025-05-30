process MULTIQC {
    tag "MultiQC report - ${file_prefix}"
    publishDir "$params.outdir/reports", mode: 'copy'
    conda "./environment.yaml"

    input:
    path html_files
    val file_prefix

    output:
    path "${file_prefix}_multiqc_report.html"

    script:
    """
    multiqc ${html_files.join(" ")} --outdir . --filename ${file_prefix}_multiqc_report.html
    """
}