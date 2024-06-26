#!/bin/bash -ue
[ ! -f  KT_1.fastq.gz ] && ln -s KT_R1.fastq.gz KT_1.fastq.gz
[ ! -f  KT_2.fastq.gz ] && ln -s KT_R2.fastq.gz KT_2.fastq.gz
trim_galore \
     \
    --cores 1 \
    --fastqc \
    --paired \
    --quality 20 \
    --length 43 \
    --adapter AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
    --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    --gzip \
    KT_1.fastq.gz \
    KT_2.fastq.gz

cat <<-END_VERSIONS > versions.yml
"TRIMGALORE":
    trimgalore: $(echo $(trim_galore --version 2>&1) | sed 's/^.*version //; s/Last.*$//')
    cutadapt: $(cutadapt --version)
END_VERSIONS
