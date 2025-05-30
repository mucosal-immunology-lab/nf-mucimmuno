process KOMPLEXITY_FILTER {
    tag "$meta.id"
    label 'process_high'
    conda "${projectDir}/environment.yaml"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}*filtered.fq.gz"), emit: reads
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        """
        [ ! -f ${prefix}.fastq.gz ] && ln -s ${reads[0]} ${prefix}.fastq.gz

        # Decompress and filter the single-end FASTQ with kz.
        zcat ${prefix}.fastq.gz | kz --filter -t 0.55 > ${prefix}_filtered.fq
        gzip -c ${prefix}_filtered.fq > ${prefix}_filtered.fq.gz
        rm ${prefix}_filtered.fq

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            kz: \$(kz --version 2>&1 | head -n 1)
        END_VERSIONS
        """
    }
    else {
        """
        [ ! -f ${prefix}_R1.fastq.gz ] && ln -s ${reads[0]} ${prefix}_R1.fastq.gz
        [ ! -f ${prefix}_R2.fastq.gz ] && ln -s ${reads[1]} ${prefix}_R2.fastq.gz

        # Decompress and filter each paired-end FASTQ with kz.
        zcat ${prefix}_R1.fastq.gz | kz --filter -t 0.55 > ${prefix}_R1_filtered.fq
        zcat ${prefix}_R2.fastq.gz | kz --filter -t 0.55 > ${prefix}_R2_filtered.fq

        # Use BBMap's repair.sh to repair the reads.
        repair.sh in1=${prefix}_R1_filtered.fq in2=${prefix}_R2_filtered.fq \
                  out1=${prefix}_R1_repaired.fq out2=${prefix}_R2_repaired.fq

        # Compress the repaired output files.
        gzip -c ${prefix}_R1_repaired.fq > ${prefix}_R1_filtered.fq.gz
        gzip -c ${prefix}_R2_repaired.fq > ${prefix}_R2_filtered.fq.gz

        # Clean up intermediate files.
        rm ${prefix}_R1_filtered.fq ${prefix}_R2_filtered.fq \
           ${prefix}_R1_repaired.fq ${prefix}_R2_repaired.fq

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            kz: \$(kz --version 2>&1 | head -n 1)
            repair.sh: \$(repair.sh 2>&1 | head -n 1)
        END_VERSIONS
        """
    }
}