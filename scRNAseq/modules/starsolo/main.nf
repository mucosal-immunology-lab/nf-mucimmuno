process STARSOLO {
    tag "STARsolo ${meta.id}"
    publishDir "$params.outdir/STARsolo/$meta.id", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta),  path('*.Solo.out')         , emit: counts
    tuple val(meta),  path('*Log.final.out')     , emit: log_final
    tuple val(meta),  path('*Log.out')           , emit: log_out
    tuple val(meta),  path('*Log.progress.out')  , emit: log_progress
    tuple val(meta),  path('*/Gene/Summary.csv') , emit: summary
    path "versions.yml"                          , emit: versions
    when:
    task.ext.when == null || task.ext.when

    script:
    def cores = 1
    if (task.cpus) {
        cores = (task.cpus as int)
    }

    // Added soft-links to original fastqs for consistent naming
    def prefix = task.ext.prefix ?: "${meta.id}"

    // Base STARsolo approach on the cell label structure (CLS) version
    switch(meta.CLS) {
        case "BD_Original":
            def whitespace_paths = ["${moduleDir}/CLS/BD_Original/BD_CLS1_original.txt",
                                    "${moduleDir}/CLS/BD_Original/BD_CLS2_original.txt",
                                    "${moduleDir}/CLS/BD_Original/BD_CLS3_original.txt"];
            def whitespace_paths_str = whitespace_paths.collect { '"' + it + '"' }.join(' ');
            solo_args = "--soloType CB_UMI_Complex " +
                "--soloCBmatchWLtype EditDist_2 " +
                "--soloCBwhitelist ${whitespace_paths_str} " +
                "--soloUMIlen 8 " +
                "--soloCBposition 0_0_0_8 0_21_0_29 0_43_0_51 " +
                "--soloUMIposition 0_52_0_59 ";
            break;
        case "BD_Enhanced_V1":
            def whitespace_paths = ["${moduleDir}/CLS/BD_Enhanced_V1/BD_CLS1_V1.txt",
                                    "${moduleDir}/CLS/BD_Enhanced_V1/BD_CLS2_V1.txt",
                                    "${moduleDir}/CLS/BD_Enhanced_V1/BD_CLS3_V1.txt"];
            def whitespace_paths_str = whitespace_paths.collect { '"' + it + '"' }.join(' ');
            solo_args = "--soloType CB_UMI_Complex " +
                "--soloCBmatchWLtype EditDist_2 " +
                "--soloCBwhitelist ${whitespace_paths_str} " +
                "--soloUMIlen 8 " +
                "--soloCBposition 2_0_2_8 2_13_2_21 3_1_3_9 " +
                "--soloUMIposition 3_10_3_17 " +
                "--soloAdapterSequence NNNNNNNNNGTGANNNNNNNNNGACA ";
            break;
        case "BD_Enhanced_V2":
            def whitespace_paths = ["${moduleDir}/CLS/BD_Enhanced_V2/BD_CLS1_V2.txt",
                                    "${moduleDir}/CLS/BD_Enhanced_V2/BD_CLS2_V2.txt",
                                    "${moduleDir}/CLS/BD_Enhanced_V2/BD_CLS3_V2.txt"];
            def whitespace_paths_str = whitespace_paths.collect { '"' + it + '"' }.join(' ');
            solo_args = "--soloType CB_UMI_Complex " +
                "--soloCBmatchWLtype EditDist_2 " +
                "--soloCBwhitelist ${whitespace_paths_str} " +
                "--soloUMIlen 8 \\" +
                "--soloCBposition 2_0_2_8 2_13_2_21 3_1_3_9 " +
                "--soloUMIposition 3_10_3_17 " +
                "--soloAdapterSequence NNNNNNNNNGTGANNNNNNNNNGACA ";
            break;
        case "10X_Chromium_V1":
            solo_args = "--soloType CB_UMI_Simple " +
                "--soloCBmatchWLtype 1MM multi Nbase pseudocounts " +
                "--soloCBwhitelist \"${moduleDir}/CLS/10X_Chromium_V1/737K-april-2014_rc.txt\" " +
                "--clipAdapterType CellRanger4 " +
                "--soloBarcodeReadLength 0 " +
                "--soloCBstart 1 " +
                "--soloCBlen 14 " +
                "--soloUMIstart 15 " +
                "--soloUMIlen 10 ";
            break;
        case "10X_Chromium_V2":
            solo_args = "--soloType CB_UMI_Simple " +
                "--soloCBmatchWLtype 1MM multi Nbase pseudocounts " +
                "--soloCBwhitelist \"${moduleDir}/CLS/10X_Chromium_V2/10X_V2_737K_august_2016.txt\" " +
                "--clipAdapterType CellRanger4 " +
                "--soloBarcodeReadLength 0 " +
                "--soloCBstart 1 " +
                "--soloCBlen 16 " +
                "--soloUMIstart 17 " +
                "--soloUMIlen 10 ";
            break;
        case "10X_Chromium_V3":
            solo_args = "--soloType CB_UMI_Simple " +
                "--soloCBmatchWLtype 1MM multi Nbase pseudocounts " +
                "--soloCBwhitelist \"${moduleDir}/CLS/10X_Chromium_V3/10X_V3_6975K_feb_2018.txt\" " +
                "--clipAdapterType CellRanger4 " +
                "--soloBarcodeReadLength 0 " +
                "--soloCBstart 1 " +
                "--soloCBlen 16 " +
                "--soloUMIstart 17 " +
                "--soloUMIlen 12 ";
            break;
        case "SeqWell":
            solo_args = "--soloType CB_UMI_Simple " +
                "--soloCBmatchWLtype 1MM multi Nbase pseudocounts " +
                "--soloCBwhitelist None " +
                "--clipAdapterType CellRanger4 " +
				"--soloStrand Unstranded " +
                "--soloBarcodeReadLength 0 " +
                "--soloCBstart 1 " +
                "--soloCBlen 12 " +
                "--soloUMIstart 13 " +
                "--soloUMIlen 8 ";
            break;
        case "DropSeq":
            solo_args = "--soloType CB_UMI_Simple " +
                "--soloCBmatchWLtype 1MM multi Nbase pseudocounts " +
                "--soloCBwhitelist None " +
                "--clipAdapterType CellRanger4 " +
				"--soloStrand Unstranded " +
                "--soloBarcodeReadLength 0 " +
                "--soloCBstart 1 " +
                "--soloCBlen 12 " +
                "--soloUMIstart 13 " +
                "--soloUMIlen 8 ";
            break;
    }

    // RUN STARsolo on the sample in either paired- or single-ended mode
    if (!meta.single_end) {
        """
        [ ! -f  ${prefix}_1.fastq.gz ] && ln -s ${reads[0]} ${prefix}_1.fastq.gz
        [ ! -f  ${prefix}_2.fastq.gz ] && ln -s ${reads[1]} ${prefix}_2.fastq.gz

        STAR \\
            --genomeDir $meta.GenomeIndex \\
            --runThreadN $cores \\
            --runMode alignReads \\
            --readFilesIn ${prefix}_2.fastq.gz ${prefix}_1.fastq.gz \\
            --readFilesCommand zcat \\
            --outFileNamePrefix ${prefix}. \\
            --soloUMIdedup $params.starsolo.soloUMIdedup \\
            --soloUMIfiltering $params.starsolo.soloUMIfiltering \\
            --soloCellFilter $params.starsolo.soloCellFilter \\
            --soloFeatures Gene GeneFull_Ex50pAS Velocyto \\
            --soloMultiMappers $params.starsolo.soloMultiMappers \\
            $solo_args \\
            --outSAMtype None
        
        if [ -d ${prefix}.Solo.out ]; then
        find ${prefix}.Solo.out \\( -name "*.tsv" -o -name "*.mtx" \\) -exec gzip {} \\;
        fi

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            star: \$(STAR --version | sed -e "s/STAR_//g")
        END_VERSIONS
        """
    } else {
        """
        [ ! -f  ${prefix}_1.fastq.gz ] && ln -s ${reads[0]} ${prefix}_1.fastq.gz
        STAR \\
            --genomeDir $meta.GenomeIndex \\
            --runThreadN $cores \\
            --runMode alignReads \\
            --readFilesIn ${prefix}_1.fastq.gz \\
            --readFilesCommand zcat \\
            --outFileNamePrefix ${prefix}. \\
            --soloType $params.starsolo.soloType \\
            --soloCBmatchWLtype $params.starsolo.soloCBmatchWLtype \\
            --soloUMIdedup $params.starsolo.soloUMIdedup \\
            --soloUMIfiltering $params.starsolo.soloUMIfiltering \\
            --soloCellFilter $params.starsolo.soloCellFilter \\
            --soloFeatures Gene GeneFull_Ex50pAS Velocyto \\
            --soloMultiMappers $params.starsolo.soloMultiMappers \\
            $solo_args \\
            --outSAMtype None
        
        if [ -d ${prefix}.Solo.out ]; then
        find ${prefix}.Solo.out \\( -name "*.tsv" -o -name "*.mtx" \\) -exec gzip {} \\;
        fi

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            star: \$(STAR --version | sed -e "s/STAR_//g")
        END_VERSIONS
        """
    }

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}.Solo.out/
    touch ${prefix}.Solo.out/Log.final.out
    touch ${prefix}.Solo.out/Log.out
    touch ${prefix}.Solo.out/Log.progress.out
    touch ${prefix}.Solo.out/Summary.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
    END_VERSIONS
    """
}
