#!/bin/bash -ue
[ ! -f  PLN_1.fastq.gz ] && ln -s PLN_R1.fastq.gz PLN_1.fastq.gz
[ ! -f  PLN_2.fastq.gz ] && ln -s PLN_R2.fastq.gz PLN_2.fastq.gz
trim_galore \
     \
    --cores 1 \
    --paired \
    --gzip \
    PLN_1.fastq.gz \
    PLN_2.fastq.gz

cat <<-END_VERSIONS > versions.yml
"TRIMGALORE":
    trimgalore: $(echo $(trim_galore --version 2>&1) | sed 's/^.*version //; s/Last.*$//')
    cutadapt: $(cutadapt --version)
END_VERSIONS
