process adapter_trimming {
    input:
    tuple val(sample_id), path(reads1), path(reads2)  // Sample ID and paired-end reads

    output:
    tuple val(sample_id), path("trimmed_R1.fastq"), path("trimmed_R2.fastq")

    errorStrategy 'ignore'  // Continue pipeline even if this process encounters errors

    script:
    """
    fastp --in1 $reads1 --in2 $reads2 --out1 trimmed_R1.fastq --out2 trimmed_R2.fastq --detect_adapter_for_pe || \\
    cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o trimmed_R1.fastq -p trimmed_R2.fastq $reads1 $reads2
    """
}
