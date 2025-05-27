process PREPARE_BRACKEN_DB {
    tag "prepare_bracken_db"
    label 'process_high'

    // Copy the final Bracken database directory out
    publishDir "$params.outdir/bracken_db", mode: 'copy'

    conda "${moduleDir}/environment.yaml"

    input:
        // Collect the Kraken2 database path
        val(kraken2_db_ch)

    output:
        // Emit the directory containing the .kmer_distrib files
        val(kraken2_db_ch), emit: bracken_db

    script:
    """
    set -euo pipefail

    echo \"Building Bracken database from Kraken2 DB: ${params.taxonomy.kraken2_db}\"

    # Run bracken-build with read length, threads, and output folder
    bracken-build \\
        -d ${kraken2_db_ch} \\
        -t ${task.cpus} \\
        -l ${params.taxonomy.kmer_length}

    echo \"Bracken database build completed successfully.\"
    """
}

process RUN_BRACKEN_CORRECTION {
    tag "$meta.id"
    label 'process_high'

    // Conda environment with bracken installed
    conda "${moduleDir}/environment.yaml"

    input:
    tuple val(meta), path(krakenReport)
    val(bracken_db_ch)

    output:
    tuple val(meta), path("${meta.id}_bracken_report.tsv"), emit: bracken_report

    script:
        // Define the sample prefix in Groovy
        def prefix = task.ext.prefix ?: "${meta.id}"
    """
    set -e

    # Run Bracken to re-estimate abundances at the desired taxonomic level
    bracken \
        -d ${bracken_db_ch} \
        -i ${krakenReport} \
        -o ${prefix}_bracken_report.tsv \
        -r ${params.taxonomy.kmer_length} \
        -l ${params.bracken.bracken_level}
    """
}

process MERGE_BRACKEN_REPORTS {
    tag "merge_bracken_reports"
    label 'process_high'
    conda "${moduleDir}/environment.yaml"

    // Publish the combined Kraken2 report into the kraken2 results folder
    publishDir "$params.outdir/bracken", mode: 'copy'

    input:
        // Collect all individual Kraken2 report files
        path bracken_reports
        val sample_ids

    output:
        // Emit the merged report
        path "combined_bracken_report.tsv",             emit: combined_report
        path "filtered_combined_bracken_report.tsv",    emit: filtered_combined_report
    
    script:
        // join into space‐separated strings for the Python call
        def namesArg = sample_ids.join(',')
        def reportsArg = bracken_reports.join(' ')
    """
    set -e

    # Run the Python combiner script from the project directory
    python ${moduleDir}/combine_bracken_reports.py \
        --files ${reportsArg} \
        --names ${namesArg} \
        -o combined_bracken_report.tsv

    # Verify the output
    [ -s combined_bracken_report.tsv ] \
        || { echo "Error: combined_bracken_report.tsv is empty" >&2; exit 1; }

    # Filter the combined report
    python3 ${moduleDir}/filter_bracken_report.py \
        -i combined_bracken_report.tsv \
        -o filtered_combined_bracken_report.tsv
    """
}