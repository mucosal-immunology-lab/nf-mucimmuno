#!/bin/bash -ue
multiqc KT_1_fastqc.html KT_2_fastqc.html SB_1_fastqc.html SB_2_fastqc.html PLN_1_fastqc.html PLN_2_fastqc.html DM_1_fastqc.html DM_2_fastqc.html KT_1_fastqc.zip KT_2_fastqc.zip SB_1_fastqc.zip SB_2_fastqc.zip PLN_1_fastqc.zip PLN_2_fastqc.zip DM_1_fastqc.zip DM_2_fastqc.zip --outdir . --filename pretrim_multiqc_report.html
