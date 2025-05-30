process PREPARE_HOST_GENOME {
    tag "PrepareHostGenome"
    publishDir "$params.outdir/hostIndex", mode: 'copy'
    conda "./environment.yaml"

    output:
        path "chm13v2.0_GRCh38_full_plus_decoy.fasta", emit: hostFasta
        path "chm13v2.0_GRCh38_full_plus_decoy.*", emit: bowtie2_index

    script:
    """
    # Check if the combined host genome file exists
    if [ ! -f chm13v2.0_GRCh38_full_plus_decoy.fasta ]; then
        echo "Downloading the CHM13 human genome..."
        wget -c https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz
        gzip -d chm13v2.0.fa.gz

        echo "Prefixing chromosome names in CHM13 genome..."
        sed 's/^>chr/>CHM13_chr/' chm13v2.0.fa > chm13v2.0_modified.fasta

        echo "Downloading the 1000 Genomes human genome..."
        wget -c http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa

        echo "Combining CHM13 and 1000 Genomes genomes..."
        cat chm13v2.0_modified.fasta GRCh38_full_analysis_set_plus_decoy_hla.fa > chm13v2.0_GRCh38_full_plus_decoy.fasta

        echo "Cleaning up intermediate files..."
        rm -rf chm13v2.0.fa chm13v2.0_modified.fasta GRCh38_full_analysis_set_plus_decoy_hla.fa

        echo "Finished preparing the human genome."
    else
        echo "Host genome file already exists. Skipping download and combination."
    fi

    echo "Indexing the host genome with bowtie2-build..."
    if ! command -v bowtie2-build &> /dev/null; then
        echo "Error: bowtie2-build is not installed or not in PATH." >&2
        exit 1
    fi
    bowtie2-build chm13v2.0_GRCh38_full_plus_decoy.fasta chm13v2.0_GRCh38_full_plus_decoy --large-index

    # Verify that the Bowtie2 index files were created successfully
    if [ ! -f "chm13v2.0_GRCh38_full_plus_decoy.1.bt2l" ]; then
        echo "Error: Bowtie2 index files were not created successfully." >&2
        exit 1
    fi
    """
}

process HOST_DECONTAMINATE {
    tag "$meta.id"
    label 'process_high'
    conda "./environment.yaml"

    // Publish decontaminated FASTQ files to $params.outdir/decontam
    publishDir "$params.outdir/decontam", mode: 'copy', pattern: "*_decontam.fq.gz"
    // Publish Bowtie2 log files to $params.outdir/decontam_log
    publishDir "$params.outdir/decontam_log", mode: 'copy', pattern: "*_bowtie2.log"

    input:
        tuple val(meta), path(reads)
        // hostIndex should be the common prefix, e.g. "/mnt/GDrive_01/.../chm13v2.0_GRCh38_full_plus_decoy"
        val(hostIndex)

    output:
        tuple val(meta), path("${meta.id}*decontam.fq.gz"), emit: reads
        path "versions.yml", emit: versions
        path "${meta.id}_bowtie2.log", emit: log

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        if (meta.single_end) {
            """
            # Run Bowtie2 for single-end reads (gzipped input, gzipped unaligned output).
            bowtie2 -p ${task.cpus} -x "${hostIndex}" -U "${reads}" --un-gz "${prefix}_decontam.fq.gz" -S /dev/null || true

            # Check that the output file exists and is non-empty.
            [ -s "${prefix}_decontam.fq.gz" ] || { echo "Error: Bowtie2 produced an empty output for single-end reads." >&2; exit 1; }

            # Copy Nextflow's generated log to a new file.
            cp .command.log "${prefix}_bowtie2.log"

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                bowtie2: \$(bowtie2 --version | head -n 1)
            END_VERSIONS
            """
        } else {
            """
            # Run Bowtie2 for paired-end reads (using --un-conc-gz to produce gzipped outputs).
            bowtie2 -p ${task.cpus} -x "${hostIndex}" -1 "${reads[0]}" -2 "${reads[1]}" --un-conc-gz "${prefix}_temp_decontam.fq.gz" -S /dev/null || true

            # Ensure that both output files exist and are not empty.
            [ -s "${prefix}_temp_decontam.fq.1.gz" ] || { echo "Error: Bowtie2 produced an empty output for R1." >&2; exit 1; }
            [ -s "${prefix}_temp_decontam.fq.2.gz" ] || { echo "Error: Bowtie2 produced an empty output for R2." >&2; exit 1; }

            # Rename the Bowtie2-produced files.
            mv "${prefix}_temp_decontam.fq.1.gz" "${prefix}_R1_decontam.fq.gz"
            mv "${prefix}_temp_decontam.fq.2.gz" "${prefix}_R2_decontam.fq.gz"

            # Copy Nextflow's log file.
            cp .command.log "${prefix}_bowtie2.log"

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                bowtie2: \$(bowtie2 --version | head -n 1)
            END_VERSIONS
            """
        }
}
