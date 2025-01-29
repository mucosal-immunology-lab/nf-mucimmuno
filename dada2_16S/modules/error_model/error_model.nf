process ERROR_MODEL {
    input:
    tuple val(run_name), path(filtered_dir)

    output:
    tuple val(run_name), path("${run_name}_filtered/error_model.rds"), emit: error_model
    path("${run_name}_filtered/error_model.pdf")
    path("${run_name}_filtered/error_model.pdf.rds"), emit: error_plots

    tag "ERROR_MODEL - ${run_name}"
    publishDir "results/${run_name}/filtering_report/", mode: 'copy'

    script:
    """
    Rscript -e '
    library(dada2)
    library(ggplot2)

    # Use number of CPUs allocated to this task for multithreading
    nc <- ${task.cpus}

    # Identify forward and reverse filtered reads
    fwd.fn <- sort(list.files("${filtered_dir}", pattern="-R1.fastq", full.names=TRUE))
    rev.fn <- sort(list.files("${filtered_dir}", pattern="-R2.fastq", full.names=TRUE))

    # Learn error models
    err <- list()
    err[[1]] <- learnErrors(fwd.fn, nbases=1e8, multithread=nc)
    err[[2]] <- learnErrors(rev.fn, nbases=1e8, multithread=nc)

    # Plot the error models
    p <- list()
    p[[1]] <- plotErrors(err[[1]], nominalQ=TRUE) + ggtitle(paste("${run_name}", "| forward reads"))
    p[[2]] <- plotErrors(err[[2]], nominalQ=TRUE) + ggtitle(paste("${run_name}", "| reverse reads"))

    # Save the error model and plots
    saveRDS(err, file.path("${filtered_dir}", "error_model.rds"))
    saveRDS(p, file.path("${filtered_dir}", "error_model.pdf.rds"))
    pdf(file.path("${filtered_dir}", "error_model.pdf"))
    invisible(lapply(p, print))
    invisible(dev.off())
    '
    """
}

process ERROR_MODEL_AGGREGATE {
    input:
    val error_plot_rds_files_list

    output:
    path "aggregated_error_plots.pdf"

    publishDir "results/", mode: 'copy'

    script:
    """
    # Write the R script into a file
    cat > aggregate_error_models.R << 'EOF'
    library(ggplot2)
    library(ggpubr)

    # List of quality RDS files
    error_model_files <- c(${error_plot_rds_files_list.collect { '"' + it + '"' }.join(', ')})

    if (length(error_model_files) == 0) {
      stop('No error_model.pdf.rds files provided')
    }

    all_plots <- lapply(error_model_files, function(file) {
      plots <- readRDS(file)
      ggarrange(plots[[1]], plots[[2]], nrow = 1)
    })

    # Save all combined plots into a single PDF, each run on a separate page
    pdf('aggregated_error_plots.pdf', width = 10, height = 5)
    invisible(lapply(all_plots, print))
    dev.off()
    EOF

    # Run the R script
    Rscript aggregate_error_models.R
    """
}