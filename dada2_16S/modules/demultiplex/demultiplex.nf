process DEMULTIPLEX {
    input:
    tuple val(run_name), path(r1_file), path(r2_file), path(index_file), path(barcode_file)

    output:
    tuple val(run_name), path("${run_name}_demultiplexed/")

    tag "DEMULTIPLEX - ${run_name}"

    script:
    """
    echo "Running DEMULTIPLEX for ${run_name}."
    mkdir -p ${run_name}_demultiplexed

    zcat ${r1_file} > ${run_name}_R1.fastq
    zcat ${r2_file} > ${run_name}_R2.fastq
    zcat ${index_file} > ${run_name}_Index.fastq

    iu-demultiplex \\
        --r1 ${run_name}_R1.fastq \\
        --r2 ${run_name}_R2.fastq \\
        -i ${run_name}_Index.fastq \\
        -x \\
        -o ${run_name}_demultiplexed \\
        -s ${barcode_file}
    """
}

process USE_EXISTING_DEMULTIPLEXED {
    when:
    demultiplexed_folder.exists()

    input:
    tuple val(run_name), path(demultiplexed_folder)

    output:
    tuple val(run_name), path("${demultiplexed_folder}")

    script:
    """
    echo "Using pre-existing demultiplexed folder for ${run_name}."
    ln -s ${demultiplexed_folder} ${run_name}_demultiplexed
    """
}