process COLLECT_EXPORT_FILES {
    tag "Collecting export files"
    publishDir "$params.outdir", mode: 'copy'

    input:
    path export_files

    output:
    path "export_files.tar.gz", emit: exported_files_tar

    script:
    """
    mkdir -p export_files

    cp -rL ${export_files.join(" ")} export_files

    find export_files -type f -print0 | parallel --null -j${task.cpus} gzip

    tar -zcvf export_files.tar.gz export_files
    """
}