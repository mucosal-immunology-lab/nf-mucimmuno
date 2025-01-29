process MERGE_PAIRED_ENDS {
    input:
    tuple val(run_name), path(filtered_dir), path(error_model_rds)
    
    output:
    tuple val(run_name), path("${run_name}_filtered/seqtab.rds"), emit: seqtab

    tag "MERGE_PAIRED_ENDS - ${run_name}"

    script:
    """
    Rscript -e '
    library(dada2)

    # Use allocated CPUs for multithreading
    nc <- ${task.cpus}

    # Load the error model passed in via input channel
    err.model <- readRDS("${error_model_rds}")

    # Identify forward and reverse filtered reads
    fwd.fn <- sort(list.files("${filtered_dir}", pattern="-R1.fastq", full.names=TRUE))
    rev.fn <- sort(list.files("${filtered_dir}", pattern="-R2.fastq", full.names=TRUE))

    # Extract sample names
    sample.names <- sapply(strsplit(basename(fwd.fn), "-R1.fastq"), `[`, 1)
    sample.names.rev <- sapply(strsplit(basename(rev.fn), "-R2.fastq"), `[`, 1)
    if(!identical(sample.names, sample.names.rev)) {
      stop("Forward and reverse files do not match.")
    }

    names(fwd.fn) <- sample.names
    names(rev.fn) <- sample.names

    # Perform merging for each sample
    merged <- vector("list", length(sample.names))
    names(merged) <- sample.names

    for(j in seq_along(sample.names)) {
      derepF <- derepFastq(fwd.fn[j])
      derepR <- derepFastq(rev.fn[j])

      asvF <- dada(derepF, err=err.model[[1]], pool=TRUE, multithread=nc)
      asvR <- dada(derepR, err=err.model[[2]], pool=TRUE, multithread=nc)

      merged[[sample.names[j]]] <- mergePairs(asvF, derepF, asvR, derepR)
    }

    # Create sequence table and save
    st <- makeSequenceTable(merged)
    saveRDS(st, file.path("${filtered_dir}", "seqtab.rds"))
    '
    """
}

process MERGE_SEQTABS {
    input:
    val seqtab_files  // A single value that is a list of file paths

    output:
    path "seqtab.rds", emit: merged_seqtab

    tag "MERGE_SEQTABS"
    publishDir "results/", mode: 'copy'

    script:
    """
    # Write the R script into a file
    cat > merge_seqtabs.R << 'EOF'
    library(dada2)
    library(dplyr)

    # File paths to sequence tables
    seqtab.fps <- c("/mnt/GDrive_01/Monash/33_GitHub/nf-mucimmuno/dada2_16S/work/17/d33ae64ccfb4dec2e14d4f5a6c2bb8/run_02_filtered/seqtab.rds", "/mnt/GDrive_01/Monash/33_GitHub/nf-mucimmuno/dada2_16S/work/e8/83a1a0c4d356843f303999592a2849/run_01_filtered/seqtab.rds")

    # Validate the presence of sequence table files
    if (length(seqtab.fps) == 0) {
      stop("No seqtab.rds files provided.")
    }

    # Function to extract run name from file path
    extract_run_name <- function(filepath) {
      folder <- dirname(filepath)  # Get the folder path
      basename(folder) %>% sub("_filtered.*\$", "", .)  # Extract run name before '_filtered'
    }

    # Process each sequence table
    seqtabs <- lapply(seqtab.fps, function(f) {
      seqtab <- readRDS(f)
      if (is.null(seqtab) || nrow(seqtab) == 0) stop(paste("Invalid or empty seqtab:", f))
      
      run_name <- extract_run_name(f)  # Extract the run name
      # Update row names to ensure uniqueness
      rownames(seqtab) <- gsub("_r.*|-r.*|_R.*|-R.*", "", rownames(seqtab))
      rownames(seqtab) <- paste0(rownames(seqtab), "_", run_name)

      seqtab
    })

    # Merge all sequence tables
    if (length(seqtabs) == 1) {
      seqtab <- seqtabs[[1]]
    } else {
      seqtab <- do.call(mergeSequenceTables, seqtabs)
    }

    # Save the merged sequence table
    saveRDS(seqtab, "seqtab.rds")

    EOF

    # Execute the R script
    Rscript merge_seqtabs.R
    """
}

process REMOVE_CHIMERAS {
    input:
    path merged_seqtab

    output:
    path "seqtab_nochim.rds", emit: seqtab_nochim
    path "seqtab_nochim.txt"
    path "length_distribution.rds"
    path "length_distribution.pdf"

    tag "REMOVE_CHIMERAS"
    publishDir "results/", mode: 'copy'

    script:
    """
    # Create an R script file
    cat > remove_chimeras.R <<EOF
    library(dada2)
    library(data.table)

    nc <- ${task.cpus}
    seqtab <- readRDS("${merged_seqtab}")

    seqtab_nochim <- removeBimeraDenovo(seqtab, method="per-sample", multithread=nc, verbose=TRUE)

    saveRDS(seqtab_nochim, "seqtab_nochim.rds")
    fwrite(as.data.frame(seqtab_nochim), "seqtab_nochim.txt", quote=FALSE, sep="\\t")

    distrib <- table(nchar(getSequences(seqtab_nochim)))
    saveRDS(distrib, "length_distribution.rds")

    pdf("length_distribution.pdf")
    plot(distrib, xlab="Read length", ylab="Number of ASVs")
    dev.off()
    EOF

    # Run the R script
    Rscript remove_chimeras.R
    """
}

process READ_TRACKING {
    input:
    val seqtab_rds
    val seqtab_nochim_rds
    val runs_info

    output:
    path "read_tracking_report.pdf"

    tag "READ_TRACKING"
    publishDir "results/", mode: 'copy'

    script:
    """
    # Write the R script into a file
    cat > read_tracking.R << 'EOF'
    library(ggplot2)
    library(ggpubr)
    library(reshape2)

    seqtab <- readRDS("${seqtab_rds}")
    seqtab_nochim <- readRDS("${seqtab_nochim_rds}")

    runs_info <- list(
      ${ runs_info.collect { run -> "c(\"${run[0]}\", \"${run[1]}\", \"${run[2]}\")" }.join(",\n  ") }
    )

    track_plots <- list()

    for (run_data in runs_info) {
      run_name <- run_data[1]
      filtered_dir <- run_data[2]
      filtering_report_path <- run_data[3]

      filtering <- readRDS(filtering_report_path)
      row.names(filtering) <- gsub("-R1.fastq", "", row.names(filtering))
      
      filtering_alt <- filtering
      rownames(filtering_alt) <- paste(rownames(filtering_alt), run_name, sep = '_')
      
      filtering <- rbind(filtering, filtering_alt)

      common_samples <- intersect(intersect(row.names(filtering), row.names(seqtab)), row.names(seqtab_nochim))
      track <- cbind(
        filtering[rownames(filtering) %in% common_samples,],
        rowSums(seqtab[rownames(seqtab) %in% common_samples, , drop=FALSE]),
        rowSums(seqtab_nochim[rownames(seqtab_nochim) %in% common_samples, , drop=FALSE])
      )
      colnames(track) <- c('Input', 'Filtered', 'Merged', 'Non-chimeric')

      # Adjust counts
      for (j in (ncol(track)-1):1) {
        for (k in (j+1):ncol(track)) {
          track[, j] <- track[, j] - track[, k]
        }
      }

      p <- ggplot(melt(as.matrix(track)), aes(x=Var1, y=value, fill=Var2)) +
           geom_col() +
           scale_fill_manual(values = c("Input" = "#d6ccc2", "Filtered" = "#2a9d8f",
                                        "Merged" = "#0077b6", "Non-chimeric" = "#5e548e"),
                             name = 'Stage') +
           theme_pubr(legend = 'right', x.text.angle = 90, base_size = 8) +
           labs(title = paste0('Reads tracking - ', run_name), x='Samples', y='Reads', fill=NULL)

      saveRDS(p, file.path(filtered_dir, 'read_tracking_report.pdf.rds'))
      pdf(file.path(filtered_dir, 'read_tracking_report.pdf'), width = 10, height = 4)
      print(p)
      dev.off()

      track_plots[[length(track_plots)+1]] <- p
    }

    # Combine all plots
    pdf("read_tracking_report.pdf", width = 10, height = 6)
    invisible(lapply(track_plots, print))
    dev.off()
  EOF

    # Run the R script
    Rscript read_tracking.R
    """
}
