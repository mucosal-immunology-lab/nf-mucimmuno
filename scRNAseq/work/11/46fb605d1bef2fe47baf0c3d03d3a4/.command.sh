#!/bin/bash -ue
printf "%s %s\n" PLN_R1.fastq.gz PLN_1.gz PLN_R2.fastq.gz PLN_2.gz | while read old_name new_name; do
    [ -f "${new_name}" ] || ln -s $old_name $new_name
done

fastqc \
     \
    --threads 2 \
    --memory 8192 \
    PLN_1.gz PLN_2.gz

cat <<-END_VERSIONS > versions.yml
"FASTQC":
    fastqc: $( fastqc --version | sed '/FastQC v/!d; s/.*v//' )
END_VERSIONS
