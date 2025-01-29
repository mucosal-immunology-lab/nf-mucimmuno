process DE_NOVO_PHYLO_TREE {
    input:
    path seqtab_nochim
    path taxonomy_species

    output:
    path "tree.rds"

    tag "DE_NOVO_PHYLO_TREE"
    publishDir "results/", mode: 'copy'

    script:
    """
    # Part 1: Write the R script for preprocessing before RAxML
    cat > preprocess_raxml.R << 'EOF'
    library(dplyr)
    library(DECIPHER)
    library(ips)
    library(phangorn)

    # Load required datasets
    seqtab_nochim <- readRDS('${seqtab_nochim}')
    taxonomy_species <- readRDS('${taxonomy_species}')

    # Extract ASV sequences
    seqs <- colnames(seqtab_nochim)
    names(seqs) <- seqs

    # Align sequences
    alignment <- AlignSeqs(DNAStringSet(seqs), anchor = NA)
    phang_align <- phyDat(as(alignment, 'matrix'), type = 'DNA')

    # Convert phyDat class to DNAbin format
    alignment_dnabin <- as.DNAbin(phang_align)

    # Prepare suitable names for RAxML
    seq_names_tidy <- data.frame(original = rownames(alignment_dnabin)) %>%
        mutate(raxml_names = 'seq') %>%
        mutate(raxml_names = gsub('\\\\.', '_', make.unique(raxml_names)))
    rownames(alignment_dnabin) <- seq_names_tidy\$raxml_names

    # Write alignment to file
    write.dna(alignment_dnabin, file = "dada2.phy", format = "interleaved")

    # Save seq_names_tidy for later correction of tip labels
    saveRDS(seq_names_tidy, "seq_names_tidy.rds")
    EOF

    # Run the preprocessing R script
    Rscript preprocess_raxml.R

    # Part 2: Call RAxML directly
    raxml_exec=\$(command -v raxmlHPC-PTHREADS)
    \$raxml_exec -T ${task.cpus} -f a -p 12345 -x 12345 -m GTRGAMMAI -N 100 -s dada2.phy -n dada2

    # Part 3: Write the R script for postprocessing after RAxML
    cat > postprocess_raxml.R << 'EOF'
    library(dplyr)
    library(ape)

    # Load seq_names_tidy
    seq_names_tidy <- readRDS("seq_names_tidy.rds")

    # Read in the resulting tree from the RAxML output file
    tree <- read.tree("RAxML_bestTree.dada2")

    # Correct the tip labels
    tip_labels <- data.frame(raxml_names = tree\$tip.label) %>%
        left_join(seq_names_tidy, by = 'raxml_names')
    tree\$tip.label <- tip_labels\$original

    # Save the corrected tree to disk
    saveRDS(tree, "tree.rds")
    EOF

    # Run the postprocessing R script
    Rscript postprocess_raxml.R
    """
}
