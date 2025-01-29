process ASSIGN_TAXONOMY {
    input:
    path seqtab_nochim

    output:
    path "taxonomy.rds", emit: taxonomy_rds
    path "taxonomy.txt", emit: taxonomy_txt

    tag "ASSIGN_TAXONOMY"
    publishDir "results/", mode: 'copy'

    script:
    """
    # Variables
    REF_URL=$params.assign_taxonomy.trainSet_link
    REF_FILE="${moduleDir}/$params.assign_taxonomy.trainSet_file"
    
    # Download the reference database if not present
    if [ ! -f \$REF_FILE ]; then
        echo "Downloading Silva reference database..."
        wget -O \$REF_FILE \$REF_URL
    else
        echo "Silva reference database already exists at \$REF_FILE. Skipping download."
    fi

    # Write the R script into a file
    cat > assign_taxonomy.R << EOF
    library(dada2)
    library(data.table)
    
    # Load the sequence table
    seqtab_nochim <- readRDS('${seqtab_nochim}')
    
    # Path to the DADA2-formatted reference database (gzipped)
    db_fp <- '\$REF_FILE'
    
    # Assign taxonomy
    taxonomy <- assignTaxonomy(seqtab_nochim, db_fp, minBoot=50, tryRC=TRUE, multithread=${task.cpus})
    
    # Save taxonomy results as RDS
    saveRDS(taxonomy, 'taxonomy.rds')
    
    # Save taxonomy results as a tab-delimited text file
    fwrite(as.data.frame(taxonomy), 'taxonomy.txt', quote=F, sep='\t')
    EOF

    # Run the Rscript
    Rscript assign_taxonomy.R
    """
}

process ASSIGN_SPECIES {
    input:
    path taxonomy

    output:
    path "taxonomy_species.rds", emit: taxonomy_species_rds
    path "taxonomy_species.txt", emit: taxonomy_species_txt

    tag "ASSIGN_SPECIES"
    publishDir "results/", mode: 'copy'

    script:
    """
    # Variables
    REF_URL=$params.assign_taxonomy.assignSpecies_link
    REF_FILE="${moduleDir}/$params.assign_taxonomy.assignSpecies_file"
    
    # Download the reference database if not present
    if [ ! -f \$REF_FILE ]; then
        echo "Downloading Silva species assignment reference database..."
        wget -O \$REF_FILE \$REF_URL
    else
        echo "Silva species assignment reference database already exists at \$REF_FILE. Skipping download."
    fi

    # Write the R script into a file
    cat > assign_species.R << EOF
    library(dada2)
    library(data.table)
    
    # Load the sequence table
    taxonomy <- readRDS('${taxonomy}')
    
    # Path to the DADA2-formatted reference database (gzipped)
    db_fp <- '\$REF_FILE'
    
    # Run the taxonomy assignment incrementally to avoid memory overloading
    chunk_size <- 10000
    taxonomy_species <- do.call(rbind, lapply(split(c(1:nrow(taxonomy)), 
                                              sort(c(1:nrow(taxonomy))%%ceiling(nrow(taxonomy)/chunk_size))),
                                              function(x){
                                                return(addSpecies(taxonomy[x, ], db_fp))
                                              }))

    # Save taxonomy results as RDS
    saveRDS(taxonomy_species, 'taxonomy_species.rds')
    
    # Save taxonomy results as a tab-delimited text file
    fwrite(as.data.frame(taxonomy_species), 'taxonomy_species.txt', quote=F, sep='\t')
    EOF

    # Run the Rscript
    Rscript assign_species.R
    """
}