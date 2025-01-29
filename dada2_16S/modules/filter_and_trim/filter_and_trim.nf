process FILTER_AND_TRIM {
    input:
    tuple val(run_name), path(demultiplexed_dir)

    output:
    tuple val(run_name), path("${run_name}_filtered/"), emit: filtered_files
    path("${run_name}_filtered/filtering_report.pdf")
    path("${run_name}_filtered/filtering_report.pdf.rds"), emit: filter_plots
    path("${run_name}_filtered/filtering_report.rds"), emit: filtering_report

    tag "FILTER_AND_TRIM - ${run_name}"
    publishDir "results/${run_name}/filtering_report/", mode: 'copy'

    script:
    """
    # Write the R script into a file
    cat > filter_and_trim.R << 'EOF'
    library(dada2)
    library(ggplot2)
    library(reshape2)
    library(dplyr)
    library(ggpubr)

    # Read in the parameters directly from Nextflow params
    truncLen <- as.numeric(unlist(strsplit("${params.filter_and_trim.truncLen}", ",")))
    maxEE <- as.numeric(unlist(strsplit("${params.filter_and_trim.maxEE}", ",")))
    trimLeft <- as.numeric(unlist(strsplit("${params.filter_and_trim.trimLeft}", ",")))
    truncQ <- as.numeric(unlist(strsplit("${params.filter_and_trim.truncQ}", ",")))
    maxN <- as.numeric(unlist(strsplit("${params.filter_and_trim.maxN}", ",")))

    # File paths for forward and reverse reads
    fwd.fn <- sort(list.files("${demultiplexed_dir}", pattern="R1.fastq", full.names=TRUE))
    rev.fn <- sort(list.files("${demultiplexed_dir}", pattern="R2.fastq", full.names=TRUE))

    # Create filtered output directory
    dir.create("${run_name}_filtered", showWarnings = FALSE)

    # Set output file paths for filtered reads
    filt.fwd <- file.path("${run_name}_filtered", basename(fwd.fn))
    filt.rev <- file.path("${run_name}_filtered", basename(rev.fn))

    # Filter and trim reads
    filterAndTrim.out <- filterAndTrim(
      fwd = fwd.fn,
      filt = filt.fwd,
      rev = rev.fn,
      filt.rev = filt.rev,
      truncLen = truncLen,
      trimLeft = trimLeft,
      maxEE = maxEE,
      truncQ = truncQ,
      maxN = maxN,
      rm.phix = TRUE,
      compress = TRUE,
      verbose = TRUE
    )

    # Save the filtering summary
    saveRDS(filterAndTrim.out, "${run_name}_filtered/filtering_report.rds")

    # Prepare data for plotting
    data <- as.data.frame(filterAndTrim.out)
    row.names(data) <- gsub("-R1.fastq", "", row.names(data))
    data <- data %>%
        mutate(reads.in = reads.in - reads.out)
    data_melted <- melt(as.matrix(data))

    # Create the plot
    p <- ggplot(data_melted, aes(x = Var1, y = value, fill = Var2)) +
         geom_col() +
         scale_fill_manual(values = c("reads.in" = "#d6ccc2", "reads.out" = "#2a9d8f")) +
         labs(title = "Filtering report - ${run_name}", x = "Samples", y = "Reads") +
         theme_pubr(legend = "right", x.text.angle = 90, base_size = 8)

    # Save the plot as both RDS and PDF
    saveRDS(p, "${run_name}_filtered/filtering_report.pdf.rds")
    pdf("${run_name}_filtered/filtering_report.pdf", width = 10, height = 4)
    print(p)
    dev.off()
    EOF

    # Run the R script
    Rscript filter_and_trim.R
    """
}

process FILTER_AND_TRIM_AGGREGATE {
    input:
    val filter_plot_rds_files_list

    output:
    path "aggregated_filtering_plots.pdf"

    publishDir "results/", mode: 'copy'

    script:
    """
    # Write the R script into a file
    cat > aggregate_filtering.R << 'EOF'
    library(ggplot2)
    library(ggpubr)

    # List of quality RDS files
    filtering_files <- c(${filter_plot_rds_files_list.collect { '"' + it + '"' }.join(', ')})

    if (length(filtering_files) == 0) {
      stop('No filtering_report.pdf.rds files provided')
    }

    all_plots <- lapply(filtering_files, function(file) {
      plots <- readRDS(file)
    })

    # Save all combined plots into a single PDF, each run on a separate page
    pdf('aggregated_filtering_plots.pdf', width = 10, height = 4)
    invisible(lapply(all_plots, print))
    dev.off()
    EOF

    # Run the R script
    Rscript aggregate_filtering.R
    """
}