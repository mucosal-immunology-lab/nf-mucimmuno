process TRIMGALORE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*{3prime,5prime,trimmed,val}*.fq.gz"), emit: reads
    tuple val(meta), path("*report.txt")                        , emit: log     , optional: true
    tuple val(meta), path("*unpaired*.fq.gz")                   , emit: unpaired, optional: true
    tuple val(meta), path("*.html")                             , emit: html    , optional: true
    tuple val(meta), path("*.zip")                              , emit: zip     , optional: true
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // Calculate number of --cores for TrimGalore based on value of task.cpus
    // See: https://github.com/FelixKrueger/TrimGalore/blob/master/Changelog.md#version-060-release-on-1-mar-2019
    // See: https://github.com/nf-core/atacseq/pull/65
    def cores = 1
    if (task.cpus) {
        cores = (task.cpus as int) - 4
        if (meta.single_end) cores = (task.cpus as int) - 3
        if (cores < 1) cores = 1
        if (cores > 8) cores = 8
    }

    // Base trimming length settings on the CLS type (and associated R1 lengths)
    def trimgalore_length = 0
    switch(meta.CLS) {
        case "BD_Original":
            trimgalore_length = 43;
            break;
        case "BD_Enhanced_V1":
            trimgalore_length = 43;
            break;
        case "BD_Enhanced_V2":
            trimgalore_length = 43;
            break;
        case "10X_Chromium_V1":
            trimgalore_length = 24;
            break;
        case "10X_Chromium_V2":
            trimgalore_length = 26;
            break;
        case "10X_Chromium_V3":
            trimgalore_length = 28;
            break;
        case "SeqWell":
            trimgalore_length = 20;
            break;
        case "DropWell":
            trimgalore_length = 20;
            break;
    }

    // Added soft-links to original fastqs for consistent naming in MultiQC
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        def args_list = args.split("\\s(?=--)").toList()
        args_list.removeAll { it.toLowerCase().contains('_r2 ') }
        """
        [ ! -f  ${prefix}.fastq.gz ] && ln -s $reads ${prefix}.fastq.gz
        trim_galore \\
            ${args_list.join(' ')} \\
            --cores $cores \\
            --fastqc \\
            --quality $params.trimgalore.quality \\
            --length ${trimgalore_length} \\
            --gzip \\
            ${prefix}.fastq.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            trimgalore: \$(echo \$(trim_galore --version 2>&1) | sed 's/^.*version //; s/Last.*\$//')
            cutadapt: \$(cutadapt --version)
        END_VERSIONS
        """
    } else {
        """
        [ ! -f  ${prefix}_1.fastq.gz ] && ln -s ${reads[0]} ${prefix}_1.fastq.gz
        [ ! -f  ${prefix}_2.fastq.gz ] && ln -s ${reads[1]} ${prefix}_2.fastq.gz
        trim_galore $args \\
            --cores $cores \\
            --fastqc \\
            --paired \\
            --quality $params.trimgalore.quality \\
            --length ${trimgalore_length} \\
            --adapter $params.trimgalore.adapter \\
            --adapter2 $params.trimgalore.adapter2 \\
            --gzip \\
            ${prefix}_1.fastq.gz \\
            ${prefix}_2.fastq.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            trimgalore: \$(echo \$(trim_galore --version 2>&1) | sed 's/^.*version //; s/Last.*\$//')
            cutadapt: \$(cutadapt --version)
        END_VERSIONS
        """
    }
}
