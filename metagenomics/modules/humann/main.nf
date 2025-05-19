process PREPARE_HUMANN_DATABASES {
    tag "prepare_humann_dbs"
    label 'process_high'

    // Publish the completed database folders to results/humann_db
    publishDir "$params.outdir/humann_db", mode: 'copy'

    conda "${moduleDir}/environment.yaml"

    output:
        // Emit the base folder name so we can pass it along as a channel
        path "humann_db", emit: humann_db

    script:
    """
    set -euo pipefail

    # Create a local directory for databases
    mkdir -p humann_db

    # Download each HUMAnN3 database into humann_db,
    # updating the HUMAnN config so downstream humann picks them up.
    humann_databases --download chocophlan full humann_db --update-config yes
    humann_databases --download uniref uniref90_diamond humann_db --update-config yes
    humann_databases --download utility_mapping full humann_db --update-config yes

    # Show where HUMAnN will look for its databases
    echo 'HUMAnN will use these database folders:'
    humann_config --print | grep database_folders
    """
}