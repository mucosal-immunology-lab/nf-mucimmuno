process DECOMPRESS_EXPORT_ARCHIVE {
    publishDir "$params.outdir", mode: 'copy'

    // Define the input directory
    input:
    path export_archive

    // Define the output
    output:
    path 'export_directory'

    script:
    """
    mkdir -p export_directory
    tar -xzvf ${export_archive} --strip-components=1 -C export_directory
    
    cd export_directory || exit

    for dir in *.Solo.out; do
        if [ -d "\$dir" ]; then
            new_name="\${dir%.Solo.out}"
            mv "\$dir" "\$new_name"
        fi
    done

    ls -d */
    """
}

process SOUPX_CORRECTION {
    tag "$meta.id"
    publishDir "$params.outdir", mode: 'copy'

    // Define the input parameters
    input:
    tuple path(sample_dir), val(data_option), val(pca_dims), val(meta)

    // Define the output
    output:
    path "${sample_dir}/${data_option}/soupX_corrected" , emit: soupX_folders

    script:
    """
    Rscript ${moduleDir}/soupX_correction.R $sample_dir $data_option $pca_dims
    """
}
