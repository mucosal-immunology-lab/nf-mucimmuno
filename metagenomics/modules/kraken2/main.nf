process PREPARE_KRAKEN2_DB {
    tag "prepare_kraken2_db"
    publishDir "$params.outdir/kraken2_db", mode: 'copy'

    conda "${moduleDir}/environment.yaml"

    output:
        // Emit the completed Kraken2 database directory.
        path "kraken2_database", emit: kraken2_db

    script:
    """
    set -e

    # Define relative directories within the work directory.
    KRAKEN2_DIR="kraken2"
    KRAKEN2_DB="kraken2_database"

    echo "Cloning the Kraken2 repository into \$KRAKEN2_DIR..."
    git clone https://github.com/DerrickWood/kraken2 "\$KRAKEN2_DIR"

    cd "\$KRAKEN2_DIR"

    echo "Installing Kraken2..."
    # Pass the current working directory as the installation directory.
    bash install_kraken2.sh "\$PWD"
    
    echo "Building utilities..."
    cd src
    make
    if [ ! -f k2mask ]; then
        echo "Error: k2mask not found in src. Please ensure the build process was successful." >&2
        exit 1
    fi
    cd ..

    # Workaround: explicitly copy k2mask into the current (top-level) working directory,
    # then add that directory to PATH.
    cp "\$PWD/src/k2mask" .
    chmod +x k2mask
    export PATH="\$(pwd):\$PATH"
    
    echo "Current PATH: \$PATH"
    if ! command -v k2mask &> /dev/null; then
        echo "Error: k2mask is not on the PATH even after copying." >&2
        exit 1
    fi
    
    # Verify that kraken2-build is available.
    KRAKEN2_BUILD="\$PWD/scripts/kraken2-build"
    if [ ! -f "\$KRAKEN2_BUILD" ]; then
        echo "Error: kraken2-build not found at \$KRAKEN2_BUILD." >&2
        exit 1
    fi

    echo "Downloading taxonomy and creating Kraken2 database at \$KRAKEN2_DB..."
    "\$KRAKEN2_BUILD" --download-taxonomy --db "\$KRAKEN2_DB"
    
    echo "Checking for BLAST installation..."
    blastn -version || { echo "Error: blastn not available." >&2; exit 1; }
    
    echo "Downloading reference libraries: archaea, bacteria, fungi, viral..."
    "\$KRAKEN2_BUILD" --download-library archaea --db "\$KRAKEN2_DB"
    "\$KRAKEN2_BUILD" --download-library bacteria --db "\$KRAKEN2_DB"
    "\$KRAKEN2_BUILD" --download-library fungi --db "\$KRAKEN2_DB"
    "\$KRAKEN2_BUILD" --download-library viral --db "\$KRAKEN2_DB"
    
    echo "Downloading human genome reference..."
    "\$KRAKEN2_BUILD" --download-library human --db "\$KRAKEN2_DB"
    
    echo "Building the Kraken2 database. This may take some time..."
    "\$KRAKEN2_BUILD" --build --db "\$KRAKEN2_DB" --kmer-len "${params.taxonomy.kmer_length}"
    
    echo "Cleaning up intermediate files..."
    "\$KRAKEN2_BUILD" --clean --db "\$KRAKEN2_DB"
    
    echo "Kraken2 database build completed successfully."
    """
}

process CLASSIFY_KRAKEN2 {
    tag "$meta.id"
    label 'process_high'
    conda "${moduleDir}/environment.yaml"

    // Publish classification outputs
    publishDir "$params.outdir/kraken2", mode: 'copy', pattern: "*.{kraken,report,versions.yml}"

    input:
        tuple val(meta), path(reads)
        val(kraken2_db)

    output:
        tuple val(meta), path("${meta.id}.kraken"), emit: kraken
        tuple val(meta), path("${meta.id}.report"), emit: report
        path("versions.yml"),                       emit: versions

    script:
        // Define the sample prefix in Groovy
        def prefix = task.ext.prefix ?: "${meta.id}"

    """
    set -e

    # Common kraken2 args
    K2_ARGS="--db ${kraken2_db} --threads ${task.cpus} --use-names --output ${prefix}.kraken --report ${prefix}.report"

    if [ "${meta.single_end}" = "true" ]; then
        echo "Classifying single-end reads for sample ${prefix}..."
        kraken2 \$K2_ARGS --gzip-compressed -U ${reads[0]}
    else
        echo "Classifying paired-end reads for sample ${prefix}..."
        kraken2 \$K2_ARGS --gzip-compressed --paired -1 ${reads[0]} -2 ${reads[1]}
    fi

    # Validate outputs
    [ -s "${prefix}.kraken" ] || { echo "Error: ${prefix}.kraken is empty" >&2; exit 1; }
    [ -s "${prefix}.report" ] || { echo "Error: ${prefix}.report is empty" >&2; exit 1; }

    # Record the version
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$(kraken2 --version 2>&1 | head -n 1)
    END_VERSIONS
    """
}

process MERGE_KRAKEN2_REPORTS {
    tag "merge_kraken2_reports"
    label 'process_high'
    conda "${moduleDir}/environment.yaml"

    // Publish the combined Kraken2 report into the kraken2 results folder
    publishDir "$params.outdir/kraken2", mode: 'copy'

    input:
        // Collect all individual Kraken2 report files
        path kraken_reports
        val sample_ids

    output:
        // Emit the merged report
        path "combined_kraken2_report.tsv", emit: combined_report
        path "filtered_combined_kraken2_report.tsv", emit: filtered_report

    script:
        // join into space‐separated strings for the Python call
        def namesArg = sample_ids.join(' ')
        def reportsArg = kraken_reports.join(' ')
    """
    set -e

    # Run the Python combiner script from the project directory
    python ${moduleDir}/combine_kreports.py \
        -r ${reportsArg} \
        --sample-names ${namesArg} \
        -o combined_kraken2_report.tsv

    # Verify the output
    [ -s combined_kraken2_report.tsv ] \
        || { echo "Error: combined_kraken2_report.tsv is empty" >&2; exit 1; }

    # Filter the combined report
    python3 ${moduleDir}/filter_combined_report.py \
        -i combined_kraken2_report.tsv \
        -o filtered_combined_kraken2_report.tsv
    """
}
