process QUALITY_CHECK {
    input:
    tuple val(run_name), path(demultiplexed_dir)

    output:
    path "${run_name}_quality_score.pdf.rds"  // RDS file with quality plots

    tag "QUALITY_CHECK - ${run_name}"
    publishDir "results/${run_name}", mode: 'copy'

    script:
    """
    # Write the R script into a file
    cat > quality_check.R << 'EOF'
    library(dada2)
    library(ggplot2)
    library(ggpubr)

    # List all forward and reverse read files separately in the demultiplexed directory
    forward_reads <- list.files('$demultiplexed_dir', pattern = 'R1.*.fastq(.gz)?', full.names = TRUE, recursive = TRUE)
    reverse_reads <- list.files('$demultiplexed_dir', pattern = 'R2.*.fastq(.gz)?', full.names = TRUE, recursive = TRUE)

    # Debug: Print the read file lists to ensure they are found correctly
    print('Forward Reads:')
    print(forward_reads)
    print('Reverse Reads:')
    print(reverse_reads)

    # Check if there are files to process
    if (length(forward_reads) == 0 || length(reverse_reads) == 0) {
      stop('No R1 or R2 FASTQ files found in the specified directory: $demultiplexed_dir')
    }

    # Generate quality profile plots for forward and reverse reads
    plots <- list(
      plotQualityProfile(forward_reads, n = 1e+06, aggregate = TRUE) + ggtitle(paste('Forward reads |', '${run_name}')),
      plotQualityProfile(reverse_reads, n = 1e+06, aggregate = TRUE) + ggtitle(paste('Reverse reads |', '${run_name}'))
    )

    # Save the RDS file
    saveRDS(plots, '${run_name}_quality_score.pdf.rds')

    # Save the plots as PDF
    pdf('${run_name}_quality_score.pdf', paper = 'a4')
    invisible(lapply(plots, print))
    dev.off()
    EOF

    # Run the R script
    Rscript quality_check.R
    """
}

process QUALITY_CHECK_AGGREGATE {
    input:
    val quality_rds_files_list

    output:
    path "aggregated_quality_profiles.pdf"

    publishDir "results/", mode: 'copy'

    script:
    """
    # Write the R script into a file
    cat > aggregate_quality.R << 'EOF'
    library(ggplot2)
    library(ggpubr)

    # List of quality RDS files
    quality_files <- c(${quality_rds_files_list.collect { '"' + it + '"' }.join(', ')})

    if (length(quality_files) == 0) {
      stop('No quality_score.pdf.rds files provided')
    }

    all_plots <- lapply(quality_files, function(file) {
      plots <- readRDS(file)
      ggarrange(plots[[1]], plots[[2]], nrow = 1)
    })

    # Save all combined plots into a single PDF, each run on a separate page
    pdf('aggregated_quality_profiles.pdf', width = 8, height = 4)
    invisible(lapply(all_plots, print))
    dev.off()
    EOF

    # Run the R script
    Rscript aggregate_quality.R
    """
}
