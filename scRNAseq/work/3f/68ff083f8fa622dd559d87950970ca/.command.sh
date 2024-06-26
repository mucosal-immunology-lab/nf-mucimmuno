#!/bin/bash -ue
[ ! -f  SB_1.fastq.gz ] && ln -s SB_R1.fastq.gz SB_1.fastq.gz
[ ! -f  SB_2.fastq.gz ] && ln -s SB_R2.fastq.gz SB_2.fastq.gz
trim_galore \
     \
    --cores 1 \
    --fastqc \
    --paired \
    --gzip \
    SB_1.fastq.gz \
    SB_2.fastq.gz

cat <<-END_VERSIONS > versions.yml
"TRIMGALORE":
    trimgalore: $(echo $(trim_galore --version 2>&1) | sed 's/^.*version //; s/Last.*$//')
    cutadapt: $(cutadapt --version)
END_VERSIONS
